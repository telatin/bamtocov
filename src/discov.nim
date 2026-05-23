import std/[algorithm, math, strutils]

type
  DiscovFormula* = enum
    dfLinear
    dfGeometric

  DiscovParams* = object
    windowLength*: int
    foldLower*: float
    foldUpper*: float
    alpha*: float
    formula*: DiscovFormula

  DiscovResult* = object
    score*: float
    propWindowsCovered*: float
    propCovWithinFoldRange*: float
    numWindows*: int

const
  DefaultDiscovParams* = DiscovParams(
    windowLength: 1000,
    foldLower: 0.5,
    foldUpper: 2.0,
    alpha: 0.5,
    formula: dfLinear
  )

proc parseDiscovFormula*(value: string): DiscovFormula =
  case value.normalize
  of "linear":
    dfLinear
  of "geometric":
    dfGeometric
  else:
    raise newException(ValueError,
        "discov formula must be one of linear, geometric; got: " & value)

proc medianSorted(values: seq[int]): float =
  let n = values.len
  if n == 0:
    return 0.0

  let mid = n div 2
  if (n and 1) == 1:
    result = values[mid].float
  else:
    result = (values[mid - 1].float + values[mid].float) / 2.0

proc validateDiscovParams*(params: DiscovParams, coverageLen: int) =
  if coverageLen <= 0:
    raise newException(ValueError, "coverage vector must not be empty")
  if params.windowLength <= 0:
    raise newException(ValueError, "discov window length must be positive")
  if params.foldLower < 0 or params.foldUpper < 0 or
      params.foldLower >= params.foldUpper:
    raise newException(ValueError, "invalid discov fold range")
  if params.alpha < 0.0 or params.alpha > 1.0:
    raise newException(ValueError, "discov alpha must be in [0, 1]")

proc computeDiscov*(coverage: seq[int], params: DiscovParams): DiscovResult =
  ## Compute the Distribution of Coverage score from per-base depths.
  validateDiscovParams(params, coverage.len)

  let
    requestedWindow = params.windowLength
    effectiveWindow = min(requestedWindow, coverage.len)

  var
    numWindows = 0
    coveredWindows = 0
    start = 0

  while start < coverage.len:
    let stop = min(start + effectiveWindow, coverage.len)
    inc(numWindows)

    var hasCoverage = false
    for i in start ..< stop:
      if coverage[i] > 0:
        hasCoverage = true
        break
      elif coverage[i] < 0:
        raise newException(ValueError, "coverage values must be nonnegative")

    if hasCoverage:
      inc(coveredWindows)
    start = stop

  if numWindows > 1:
    let remainder = coverage.len mod effectiveWindow
    if remainder > 0 and remainder.float < 0.1 * requestedWindow.float:
      dec(numWindows)
      let finalStart = coverage.len - remainder
      var finalHadCoverage = false
      for i in finalStart ..< coverage.len:
        if coverage[i] > 0:
          finalHadCoverage = true
          break
      if finalHadCoverage:
        dec(coveredWindows)

  let propWindowsCovered = coveredWindows.float / numWindows.float

  var nonzero = newSeq[int]()
  for depth in coverage:
    if depth > 0:
      nonzero.add(depth)
    elif depth < 0:
      raise newException(ValueError, "coverage values must be nonnegative")

  var propCovWithinFoldRange = 0.0
  if nonzero.len > 0:
    algorithm.sort(nonzero)
    let
      med = medianSorted(nonzero)
      lowerBound = params.foldLower * med
      upperBound = params.foldUpper * med

    var withinRange = 0
    for depth in nonzero:
      let value = depth.float
      if value >= lowerBound and value <= upperBound:
        inc(withinRange)

    propCovWithinFoldRange = withinRange.float / nonzero.len.float

  let score =
    case params.formula
    of dfLinear:
      params.alpha * propWindowsCovered +
          (1.0 - params.alpha) * propCovWithinFoldRange
    of dfGeometric:
      pow(propWindowsCovered, params.alpha) *
          pow(propCovWithinFoldRange, 1.0 - params.alpha)

  DiscovResult(
    score: score,
    propWindowsCovered: propWindowsCovered,
    propCovWithinFoldRange: propCovWithinFoldRange,
    numWindows: numWindows
  )
