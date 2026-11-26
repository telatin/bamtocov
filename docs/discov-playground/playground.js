"use strict";

function summarizeWindows(depths, windowLength) {
  const width = Math.min(windowLength, depths.length);
  const windows = [];
  let coveredWindows = 0;
  let ignoredTail = 0;

  for (let start = 0; start < depths.length; start += width) {
    const end = Math.min(start + width, depths.length);
    const ignored = start > 0 && end - start < 0.1 * width;
    const mean = depths.subarray(start, end).reduce((sum, depth) => sum + depth, 0) / (end - start);
    const covered = mean > 0;
    windows.push({ start, end, ignored, covered, mean });
    if (ignored) ignoredTail = end - start;
    else if (covered) coveredWindows++;
  }

  const totalWindows = windows.length - Number(ignoredTail > 0);
  return { windows, width, coveredWindows, totalWindows, ignoredTail, spread: coveredWindows / totalWindows };
}

function combineSpreadEvenness(spread, evenness, alpha, formula) {
  return formula === "geometric"
    ? spread ** alpha * evenness ** (1 - alpha)
    : alpha * spread + (1 - alpha) * evenness;
}

// Mirrors computeDiscov in src/discov.nim for nonnegative per-base depths.
function calculateDiscov(depths, { windowLength, lower, upper, alpha, formula }) {
  const summary = summarizeWindows(depths, windowLength);
  const nonzero = depths.filter(depth => depth > 0).sort();
  const coveredBases = nonzero.length;
  const mid = Math.floor(coveredBases / 2);
  const median = coveredBases ? (nonzero[mid] + nonzero[Math.floor((coveredBases - 1) / 2)]) / 2 : 0;
  const lowerDepth = lower * median;
  const upperDepth = upper * median;
  const withinRange = nonzero.reduce((count, depth) => count + Number(depth >= lowerDepth && depth <= upperDepth), 0);
  const evenness = coveredBases ? withinRange / coveredBases : 0;
  const score = combineSpreadEvenness(summary.spread, evenness, alpha, formula);

  return { ...summary, coveredBases, median, lowerDepth, upperDepth, withinRange, evenness, score };
}

function calculateDiscov2(depths, readCount, { windowLength, alpha, formula }) {
  const summary = summarizeWindows(depths, windowLength);
  let alignedBases = 0;
  let coveredBases = 0;
  for (const depth of depths) {
    alignedBases += depth;
    if (depth > 0) coveredBases++;
  }
  const detection = coveredBases / depths.length;
  const bonus = detection ** 2;
  const expected = -Math.expm1(-(alignedBases + readCount * summary.width) / depths.length);
  const normalized = expected > 0 ? Math.min(1, summary.spread / expected) : 0;

  let cv = null;
  if (summary.coveredWindows >= 3) {
    const means = summary.windows.filter(window => window.covered && !window.ignored).map(window => window.mean);
    const mean = means.reduce((sum, value) => sum + value, 0) / means.length;
    const variance = means.reduce((sum, value) => sum + (value - mean) ** 2, 0) / means.length;
    cv = Math.sqrt(variance) / mean;
  }
  const evenness = cv === null ? 0 : 1 / (1 + cv);
  const combined = combineSpreadEvenness(normalized, evenness, alpha, formula);
  return { ...summary, alignedBases, coveredBases, detection, bonus, expected: expected || 0,
    normalized, cv, evenness, combined, score: Math.min(1, combined + bonus) };
}

const referenceLength = 10000;
const get = id => document.getElementById(id);
const track = get("alignments");
const reference = get("reference");
const canvas = get("depth");
const parameters = get("parameters");
const parameters2 = get("parameters2");
const windowInput = get("window");
const lowerInput = get("lower");
const upperInput = get("upper");
const alphaInput = get("alpha");
const readLength = get("read-length");
const readStart = get("read-start");
const addForm = get("add-form");
const reads = new Map();
const format = new Intl.NumberFormat("en", { maximumFractionDigits: 6 });
const number = value => format.format(value);
let lastId = 0;
let drag = null;
let ignoreClick = false;
let depths, calculation, peak;

function createRead(start, length) {
  const button = document.createElement("button");
  const id = String(++lastId);
  button.type = "button";
  button.className = "read";
  button.dataset.readId = id;
  reads.set(id, { id, start, length, button });
  track.append(button);
}

function removeRead(read, focus = false) {
  reads.delete(read.id);
  read.button.remove();
  refresh();
  if (focus) (reads.values().next().value?.button ?? addForm.querySelector("button")).focus();
}

function loadPreset(starts) {
  reads.clear();
  track.replaceChildren();
  for (const start of starts) createRead(start, 500);
  refresh();
}

function renderReads() {
  const rowEnds = [];
  const ordered = [...reads.values()].sort((a, b) => a.start - b.start);
  for (const read of ordered) {
    const end = read.start + read.length;
    let row = rowEnds.findIndex(rowEnd => rowEnd <= read.start);
    if (row === -1) row = rowEnds.length;
    rowEnds[row] = end;
    read.button.style.left = `${100 * read.start / referenceLength}%`;
    read.button.style.width = `${100 * read.length / referenceLength}%`;
    read.button.style.bottom = `${4 + row * 16}px`;
    const label = `Read ${read.id}: ${number(read.start)} to ${number(end)} bp (end excluded)`;
    read.button.setAttribute("aria-label", label);
    read.button.title = label;
  }
  track.style.height = `${8 + 16 * rowEnds.length}px`;
  get("read-count").textContent = `${reads.size} ${reads.size === 1 ? "read" : "reads"}`;
}

function depthVector() {
  const changes = new Int32Array(referenceLength + 1);
  for (const { start, length } of reads.values()) {
    changes[start]++;
    changes[start + length]--;
  }
  let depth = 0;
  return changes.subarray(0, referenceLength).map(change => depth += change);
}

function drawDepth() {
  const ratio = window.devicePixelRatio || 1;
  const width = Math.round(canvas.clientWidth * ratio);
  const height = Math.round(canvas.clientHeight * ratio);
  if (!width || !height) return;
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext("2d");
  const style = getComputedStyle(canvas);
  const color = name => style.getPropertyValue(`--${name}`).trim();
  ctx.fillStyle = color("uncovered");
  ctx.fillRect(0, 0, width, height);

  for (const item of calculation?.windows ?? []) {
    const left = Math.round(item.start * width / referenceLength);
    const right = Math.round(item.end * width / referenceLength);
    if (right <= left) continue;
    ctx.fillStyle = color(item.ignored ? "ignored" : item.covered ? "covered" : "uncovered");
    ctx.fillRect(left, 0, right - left, height);
    if (left > 0 && right - left >= 8 * ratio) {
      ctx.fillStyle = color("line");
      ctx.fillRect(left, 0, 1, height);
    }
  }

  if (!peak) return;
  ctx.fillStyle = color("depth");
  for (let x = 0; x < width; x++) {
    const from = Math.floor(x * referenceLength / width);
    const to = Math.ceil((x + 1) * referenceLength / width);
    let binPeak = 0;
    for (let i = from; i < to; i++) binPeak = Math.max(binPeak, depths[i]);
    const bar = Math.round(height * binPeak / peak);
    ctx.fillRect(x, height - bar, 1, bar);
  }
}

function refresh() {
  renderReads();
  depths = depthVector();
  peak = depths.reduce((max, depth) => Math.max(max, depth), 0);
  const lower = lowerInput.valueAsNumber;
  const upper = upperInput.valueAsNumber;
  const alpha = alphaInput.valueAsNumber;
  const formula = get("formula").value;
  upperInput.setCustomValidity(upper <= lower ? "The upper fold must be greater than the lower fold." : "");
  const valid = [windowInput, lowerInput, upperInput, alphaInput].every(input => input.validity.valid);
  calculation = valid ? calculateDiscov(depths, { windowLength: windowInput.valueAsNumber, lower, upper, alpha, formula }) : null;
  get("parameter-error").hidden = valid;
  get("results").hidden = !valid;
  get("depth-detail").hidden = !valid;
  get("tail-detail").hidden = !calculation?.ignoredTail;
  get("max-depth").textContent = `Peak: ${number(peak)}×`;

  if (calculation) {
    const { spread, evenness, score, coveredWindows, totalWindows, withinRange,
      coveredBases, median, lowerDepth, upperDepth, ignoredTail } = calculation;
    get("spread-value").textContent = number(spread);
    get("spread-detail").textContent = `${number(coveredWindows)} / ${number(totalWindows)} covered windows`;
    get("evenness-value").textContent = number(evenness);
    get("evenness-detail").textContent = coveredBases
      ? `${number(withinRange)} / ${number(coveredBases)} covered bases within range`
      : "No covered bases";
    get("score").textContent = number(score);
    get("equation").textContent = formula === "geometric"
      ? `S^${number(alpha)} × E^${number(1 - alpha)}`
      : `${number(alpha)} × S + ${number(1 - alpha)} × E`;
    get("depth-detail").textContent = coveredBases
      ? `Median nonzero depth: ${number(median)}×. Evenness range: ${number(lowerDepth)}× to ${number(upperDepth)}×, inclusive.`
      : "The evenness score is 0 when all depths are zero.";
    get("tail-detail").textContent = `The final ${number(ignoredTail)} bp are ignored for spread because the last window is shorter than 10% of the window size. These bases still count toward evenness.`;
  }
  drawDepth();
  renderDiscov2();
}

function renderDiscov2() {
  const percentageMode = get("window-mode2").value === "percentage";
  const window = get("window2");
  const percentage = get("percentage2");
  const minimum = get("minimum2");
  const weight = get("alpha2");
  for (const [input, label, active] of [
    [window, "fixed-options2", !percentageMode],
    [percentage, "percentage-options2", percentageMode],
    [minimum, "minimum-options2", percentageMode]
  ]) {
    input.disabled = !active;
    get(label).hidden = !active;
  }
  const valid = (percentageMode ? [percentage, minimum, weight] : [window, weight])
    .every(input => input.validity.valid);
  get("parameter-error2").hidden = valid;
  get("results2").hidden = !valid;
  get("inputs2").hidden = !valid;
  get("tail-detail2").hidden = true;
  if (!valid) return;

  const windowLength = percentageMode
    ? Math.max(minimum.valueAsNumber, Math.floor(depths.length * percentage.valueAsNumber / 100))
    : window.valueAsNumber;
  const alpha = weight.valueAsNumber;
  const formula = get("formula2").value;
  const result = calculateDiscov2(depths, reads.size, { windowLength, alpha, formula });
  get("observed2").textContent = number(result.spread);
  get("windows2").textContent = `${number(result.coveredWindows)} / ${number(result.totalWindows)} covered windows`;
  get("expected2").textContent = number(result.expected);
  get("normalized2").textContent = number(result.normalized);
  get("evenness2").textContent = number(result.evenness);
  get("cv2").textContent = result.cv === null
    ? "CV unavailable with fewer than 3 covered windows"
    : `CV = ${number(result.cv)}. E = 1 / (1 + CV)`;
  get("detection2").textContent = number(result.detection);
  get("covered-bases2").textContent = `${number(result.coveredBases)} / ${number(depths.length)} bases covered`;
  get("bonus2").textContent = number(result.bonus);
  get("combined2").textContent = number(result.combined);
  get("equation2").textContent = formula === "geometric"
    ? `S_norm^${number(alpha)} × E^${number(1 - alpha)}`
    : `${number(alpha)} × S_norm + ${number(1 - alpha)} × E`;
  get("score2").textContent = number(result.score);
  get("inputs2").textContent = `G = ${number(depths.length)} bp, N = ${number(reads.size)} reads, B = ${number(result.alignedBases)} aligned bases, w = ${number(result.width)} bp.`;
  if (result.ignoredTail) {
    get("tail-detail2").hidden = false;
    get("tail-detail2").textContent = `The final ${number(result.ignoredTail)} bp are excluded from the window statistics. They still count toward aligned bases and detection.`;
  }
}

function readAt(target) {
  return reads.get(target.closest(".read")?.dataset.readId);
}

function baseAt(clientX) {
  const bounds = track.getBoundingClientRect();
  return (clientX - bounds.left) * referenceLength / bounds.width;
}

function clampStart(start, length) {
  return Math.min(referenceLength - length, Math.max(0, Math.round(start)));
}

function addReadAtPointer(event) {
  if (ignoreClick || readAt(event.target) || !readLength.reportValidity()) return;
  const length = readLength.valueAsNumber;
  createRead(clampStart(baseAt(event.clientX), length), length);
  refresh();
}

track.addEventListener("click", addReadAtPointer);
reference.addEventListener("click", addReadAtPointer);
reference.addEventListener("pointerdown", () => { ignoreClick = false; });

track.addEventListener("dblclick", event => {
  const read = readAt(event.target);
  if (read && !ignoreClick) removeRead(read);
});

track.addEventListener("pointerdown", event => {
  if (!event.isPrimary || event.button !== 0) return;
  ignoreClick = false;
  const read = readAt(event.target);
  if (!read) return;
  drag = { read, id: event.pointerId, x: event.clientX, start: read.start };
  read.button.setPointerCapture(event.pointerId);
});

track.addEventListener("pointermove", event => {
  if (drag?.id !== event.pointerId) return;
  if (!ignoreClick && Math.abs(event.clientX - drag.x) < 3) return;
  ignoreClick = true;
  const { read, start, x } = drag;
  read.button.classList.add("dragging");
  const next = clampStart(start + baseAt(event.clientX) - baseAt(x), read.length);
  if (read.start === next) return;
  read.start = next;
  refresh();
});

function endDrag(event) {
  if (drag?.id !== event.pointerId) return;
  drag.read.button.classList.remove("dragging");
  drag = null;
}

for (const event of ["pointerup", "pointercancel", "lostpointercapture"]) {
  track.addEventListener(event, endDrag);
}

track.addEventListener("keydown", event => {
  const read = readAt(event.target);
  if (!read) return;
  if (event.key === "Delete" || event.key === "Backspace") {
    event.preventDefault();
    removeRead(read, true);
  } else if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
    event.preventDefault();
    const direction = event.key === "ArrowLeft" ? -1 : 1;
    read.start = clampStart(read.start + direction * (event.shiftKey ? 100 : 10), read.length);
    refresh();
  }
});

readLength.addEventListener("input", () => {
  if (!readLength.validity.valid) return;
  readStart.max = referenceLength - readLength.valueAsNumber;
  if (readStart.valueAsNumber > Number(readStart.max)) readStart.value = readStart.max;
});

addForm.addEventListener("submit", event => {
  event.preventDefault();
  createRead(readStart.valueAsNumber, readLength.valueAsNumber);
  refresh();
});

parameters.addEventListener("input", refresh);
parameters.addEventListener("submit", event => event.preventDefault());
parameters.addEventListener("reset", () => requestAnimationFrame(refresh));
parameters2.addEventListener("input", renderDiscov2);
parameters2.addEventListener("submit", event => event.preventDefault());
parameters2.addEventListener("reset", () => requestAnimationFrame(renderDiscov2));
get("spread").addEventListener("click", () => loadPreset(Array.from({ length: 10 }, (_, i) => i * 1000 + 250)));
get("pileup").addEventListener("click", () => loadPreset(Array.from({ length: 10 }, () => 4250)));
get("clear").addEventListener("click", () => loadPreset([]));

loadPreset([250, 1250, 1650, 3300, 3600, 3800, 4250, 5250, 6500, 8300]);
new ResizeObserver(drawDepth).observe(canvas);
new MutationObserver(drawDepth).observe(document.documentElement, {
  attributes: true,
  attributeFilter: ["data-theme"]
});
