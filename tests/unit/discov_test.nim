import std/[math, unittest]

import ../../src/discov

proc closeEnough(a, b: float): bool =
  abs(a - b) < 0.000001

suite "DisCov":
  test "all-zero coverage returns zero components and score":
    let result = computeDiscov(@[0, 0, 0, 0], DiscovParams(
      windowLength: 2,
      foldLower: 0.5,
      foldUpper: 2.0,
      alpha: 0.5,
      formula: dfLinear
    ))

    check result.numWindows == 2
    check closeEnough(result.propWindowsCovered, 0.0)
    check closeEnough(result.propCovWithinFoldRange, 0.0)
    check closeEnough(result.score, 0.0)

  test "uniform covered vector scores one":
    let result = computeDiscov(@[3, 3, 3, 3], DiscovParams(
      windowLength: 2,
      foldLower: 0.5,
      foldUpper: 2.0,
      alpha: 0.5,
      formula: dfLinear
    ))

    check result.numWindows == 2
    check closeEnough(result.propWindowsCovered, 1.0)
    check closeEnough(result.propCovWithinFoldRange, 1.0)
    check closeEnough(result.score, 1.0)

  test "linear formula combines spread and evenness":
    let result = computeDiscov(@[5, 0, 0, 0], DiscovParams(
      windowLength: 2,
      foldLower: 0.5,
      foldUpper: 2.0,
      alpha: 0.5,
      formula: dfLinear
    ))

    check result.numWindows == 2
    check closeEnough(result.propWindowsCovered, 0.5)
    check closeEnough(result.propCovWithinFoldRange, 1.0)
    check closeEnough(result.score, 0.75)

  test "short vectors use a single whole-reference window":
    let result = computeDiscov(@[1, 0], DiscovParams(
      windowLength: 1000,
      foldLower: 0.5,
      foldUpper: 2.0,
      alpha: 0.5,
      formula: dfLinear
    ))

    check result.numWindows == 1
    check closeEnough(result.propWindowsCovered, 1.0)
    check closeEnough(result.score, 1.0)

  test "tiny trailing windows are dropped":
    var coverage = newSeq[int](41)
    coverage[40] = 10
    let result = computeDiscov(coverage, DiscovParams(
      windowLength: 20,
      foldLower: 0.5,
      foldUpper: 2.0,
      alpha: 0.5,
      formula: dfLinear
    ))

    check result.numWindows == 2
    check closeEnough(result.propWindowsCovered, 0.0)
    check closeEnough(result.score, 0.5)

  test "larger trailing windows are kept":
    var coverage = newSeq[int](43)
    coverage[40] = 10
    let result = computeDiscov(coverage, DiscovParams(
      windowLength: 20,
      foldLower: 0.5,
      foldUpper: 2.0,
      alpha: 0.5,
      formula: dfLinear
    ))

    check result.numWindows == 3
    check closeEnough(result.propWindowsCovered, 1.0 / 3.0)
    check closeEnough(result.score, 2.0 / 3.0)

  test "geometric formula is stricter than linear when spread is partial":
    let result = computeDiscov(@[5, 0, 0, 0], DiscovParams(
      windowLength: 2,
      foldLower: 0.5,
      foldUpper: 2.0,
      alpha: 0.5,
      formula: dfGeometric
    ))

    check closeEnough(result.score, sqrt(0.5))

  test "formula parser rejects unknown formulas":
    expect ValueError:
      discard parseDiscovFormula("median")

