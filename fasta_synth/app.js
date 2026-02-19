const rangeOutputs = [
  "noteLength",
  "filterCutoff",
  "resonance",
  "attack",
  "decay",
  "sustain",
  "release",
  "lfoRate",
  "lfoDepth",
  "delayTime",
  "delayFeedback",
  "delayMix"
];

rangeOutputs.forEach((id) => {
  const input = document.getElementById(id);
  const output = document.getElementById(`${id}Value`);
  if (!input || !output) return;
  const formatValue = (val) => {
    const num = Number(val);
    if (!Number.isFinite(num)) return val;
    if (input.step && input.step.includes(".")) {
      return Number(num.toFixed(3)).toString();
    }
    return Math.round(num).toString();
  };
  const updateOutput = () => {
    output.textContent = formatValue(input.value);
  };
  input.addEventListener("input", updateOutput);
  updateOutput();
});

const statusEl = document.getElementById("status");
const playBtn = document.getElementById("playBtn");
const stopBtn = document.getElementById("stopBtn");
const fastaInput = document.getElementById("fastaInput");
const maxNotesInput = document.getElementById("maxNotes");
const keySelect = document.getElementById("keySelect");
const octaveMinSelect = document.getElementById("octaveMin");
const octaveMaxSelect = document.getElementById("octaveMax");
const waveformSelect = document.getElementById("waveformSelect");
const lfoRateInput = document.getElementById("lfoRate");
const lfoDepthInput = document.getElementById("lfoDepth");
const delayTimeInput = document.getElementById("delayTime");
const delayFeedbackInput = document.getElementById("delayFeedback");
const delayMixInput = document.getElementById("delayMix");
const sequenceStartInput = document.getElementById("sequenceStart");
const sequenceEndInput = document.getElementById("sequenceEnd");
const reloadDefaultBtn = document.getElementById("reloadDefaultBtn");
const fastaFileInput = document.getElementById("fastaFile");
const fastaStatus = document.getElementById("fastaStatus");
const mobileReminder = document.getElementById("mobileReminder");
const dismissReminderBtn = document.getElementById("dismissReminderBtn");
const defaultPlayLabel = playBtn?.textContent?.trim() || "Play Sequence";
let reminderActive = false;
const helixContainer = document.getElementById("helixVisualizer");
const controls = {
  noteLength: document.getElementById("noteLength"),
  filterCutoff: document.getElementById("filterCutoff"),
  resonance: document.getElementById("resonance"),
  attack: document.getElementById("attack"),
  decay: document.getElementById("decay"),
  sustain: document.getElementById("sustain"),
  release: document.getElementById("release"),
  lfoRate: lfoRateInput,
  lfoDepth: lfoDepthInput,
  delayTime: delayTimeInput,
  delayFeedback: delayFeedbackInput,
  delayMix: delayMixInput,
};

const MINOR_PENT_INTERVALS = [0, 3, 5, 7, 10];
const NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"];
const KEY_ROOTS = {
  C: 48,
  "C#": 49,
  D: 50,
  "D#": 51,
  E: 52,
  F: 53,
  "F#": 54,
  G: 55,
  "G#": 56,
  A: 57,
  "A#": 58,
  B: 59,
};

const baseMap = { A: 0, C: 1, G: 2, T: 3, U: 3, N: 4, R: 5, Y: 6, M: 7, K: 8, S: 9 };

const BASE_COLORS = {
  A: "#f97316",
  C: "#38bdf8",
  G: "#a3e635",
  T: "#facc15",
  U: "#facc15",
  N: "#e2e8f0",
  R: "#fb7185",
  Y: "#4ade80",
  M: "#c084fc",
  K: "#38bdf8",
  S: "#f472b6",
};
const COMPLEMENT_BASES = {
  A: "T",
  T: "A",
  U: "A",
  C: "G",
  G: "C",
  R: "Y",
  Y: "R",
  M: "K",
  K: "M",
  S: "S",
};
const HELIX_SLOT_COUNT = 8;
let helixSlots = [];
let helixQueue = [];

const isMobileDevice = () => /Mobi|Android|iPhone|iPad/i.test(navigator.userAgent || "");

const parsePositiveInt = (value, fallback) => {
  const num = Number(value);
  if (Number.isFinite(num) && num > 0) {
    return Math.floor(num);
  }
  return fallback;
};

const setPlayButtonState = (disabled, label) => {
  if (!playBtn) return;
  playBtn.disabled = !!disabled;
  if (label) {
    playBtn.textContent = label;
  } else {
    playBtn.textContent = defaultPlayLabel;
  }
};

const showMobileReminder = () => {
  if (!mobileReminder) return;
  mobileReminder.classList.remove("hidden");
  reminderActive = true;
  setPlayButtonState(true, "Dismiss reminder to play");
  reportStatus("Flip silent mode off, then dismiss the reminder to start playback.");
};

const dismissMobileReminder = () => {
  if (!mobileReminder) return;
  mobileReminder.classList.add("hidden");
  reminderActive = false;
  setPlayButtonState(false);
  reportStatus("Silent-mode reminder dismissed. Tap Play when you're ready.");
};

const initHelix = () => {
  if (!helixContainer) return;
  helixContainer.innerHTML = "";
  helixSlots = [];
  helixQueue = [];
  for (let i = 0; i < HELIX_SLOT_COUNT; i += 1) {
    const rung = document.createElement("div");
    rung.className = "helix-rung";
    rung.setAttribute("aria-hidden", "true");

    const left = document.createElement("span");
    left.className = "helix-strand left";
    left.textContent = "·";

    const connector = document.createElement("div");
    connector.className = "helix-connector";
    const noteSpan = document.createElement("span");
    noteSpan.className = "helix-note";
    noteSpan.textContent = "--";
    connector.appendChild(noteSpan);

    const right = document.createElement("span");
    right.className = "helix-strand right";
    right.textContent = "·";

    rung.appendChild(left);
    rung.appendChild(connector);
    rung.appendChild(right);

    helixContainer.appendChild(rung);
    helixSlots.push({ root: rung, left, right, note: noteSpan });
  }
};

const resetHelix = () => {
  helixQueue = [];
  helixSlots.forEach((slot) => {
    slot.root.classList.remove("active");
    slot.root.style.removeProperty("--base-color");
    slot.root.style.removeProperty("--compl-color");
    slot.left.textContent = "·";
    slot.right.textContent = "·";
    slot.note.textContent = "--";
  });
};

const renderHelix = () => {
  helixSlots.forEach((slot, idx) => {
    const data = helixQueue[idx];
    if (data) {
      slot.root.classList.add("active");
      slot.root.style.setProperty("--base-color", data.color);
      slot.root.style.setProperty("--compl-color", data.complementColor);
      slot.left.textContent = data.base;
      slot.right.textContent = data.complement;
      slot.note.textContent = data.note;
    } else {
      slot.root.classList.remove("active");
      slot.root.style.removeProperty("--base-color");
      slot.root.style.removeProperty("--compl-color");
      slot.left.textContent = "·";
      slot.right.textContent = "·";
      slot.note.textContent = "--";
    }
  });
};

const getComplementBase = (base) => COMPLEMENT_BASES[base] || base || "N";

const pushHelixEvent = (base, midi) => {
  if (!helixSlots.length) return;
  const normalizedBase = (base || "N").toUpperCase();
  const complement = getComplementBase(normalizedBase);
  const color = BASE_COLORS[normalizedBase] || "#f8fafc";
  const complementColor = BASE_COLORS[complement] || "#f1f5f9";
  helixQueue.unshift({
    base: normalizedBase,
    complement,
    midi,
    note: formatMidiNote(midi),
    color,
    complementColor,
  });
  if (helixQueue.length > HELIX_SLOT_COUNT) {
    helixQueue.pop();
  }
  renderHelix();
};

const getRootMidi = () => {
  const key = keySelect?.value || "A";
  return KEY_ROOTS[key] ?? KEY_ROOTS.A;
};

const getOctaveLimits = () => {
  const minRaw = Number(octaveMinSelect?.value ?? 48);
  const maxRaw = Number(octaveMaxSelect?.value ?? 72);
  return {
    min: Math.min(minRaw, maxRaw),
    max: Math.max(minRaw, maxRaw),
  };
};

const getOctaveRangeLabel = () => {
  const { min, max } = getOctaveLimits();
  return `${formatMidiNote(min)}–${formatMidiNote(max)}`;
};

const formatMidiNote = (midi) => {
  const octave = Math.floor(midi / 12) - 1;
  const name = NOTE_NAMES[midi % 12] ?? "C";
  return `${name}${octave}`;
};

const clampToOctaveWindow = (midi) => {
  const { min, max } = getOctaveLimits();
  if (Number.isNaN(min) || Number.isNaN(max)) return midi;
  if (min >= max) return min;
  let note = midi;
  const guard = 20;
  let iter = 0;
  while (note > max && iter < guard) {
    note -= 12;
    iter += 1;
  }
  while (note < min && iter < guard) {
    note += 12;
    iter += 1;
  }
  if (note < min || note > max) {
    const midpoint = min + (max - min) / 2;
    note = Math.round(midpoint);
  }
  return note;
};

const DEFAULT_FASTA_PATH = "XIST.fasta";

const normalizeSequenceWindowInputs = () => {
  const startVal = Math.max(1, parsePositiveInt(sequenceStartInput?.value, 1));
  let endVal = parsePositiveInt(sequenceEndInput?.value, startVal);
  if (!Number.isFinite(endVal) || endVal < startVal) {
    endVal = startVal;
  }
  if (sequenceStartInput) {
    sequenceStartInput.value = String(startVal);
  }
  if (sequenceEndInput) {
    sequenceEndInput.value = String(endVal);
  }
  return { start: startVal, end: endVal };
};

const getSequenceWindow = (sequence) => {
  if (!sequence.length) {
    return { sequence: "", start: 0, end: 0 };
  }
  const { start: normalizedStart, end: normalizedEnd } = normalizeSequenceWindowInputs();
  const fallbackEnd = Math.min(sequence.length, normalizedEnd || normalizedStart);
  const parsedStart = normalizedStart;
  let parsedEnd = normalizedEnd;
  if (!Number.isFinite(parsedEnd) || parsedEnd <= 0) {
    parsedEnd = fallbackEnd;
  }
  const startVal = Math.max(1, parsedStart);
  let endVal = Math.max(startVal, parsedEnd);
  if (!Number.isFinite(endVal) || endVal <= 0) {
    endVal = Math.min(sequence.length, startVal + 999);
  }
  const startIdx = Math.min(sequence.length - 1, startVal - 1);
  const endIdx = Math.min(sequence.length, Math.max(startIdx + 1, endVal));
  return {
    sequence: sequence.slice(startIdx, endIdx),
    start: startIdx + 1,
    end: endIdx,
  };
};

const setFastaStatus = (message, color = "#94a3b8") => {
  if (!fastaStatus) return;
  fastaStatus.textContent = message;
  fastaStatus.style.color = color;
};

const loadDefaultFasta = async (announce = true) => {
  if (!fastaInput) return;
  setFastaStatus("Loading XIST.fasta...");
  try {
    const response = await fetch(DEFAULT_FASTA_PATH, { cache: "no-store" });
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}`);
    }
    const text = await response.text();
    if (!text.trim()) {
      throw new Error("XIST.fasta was empty");
    }
    fastaInput.value = text.trim();
    setFastaStatus("Loaded XIST.fasta reference.", "#22c55e");
    if (announce) {
      reportStatus("Loaded default XIST.fasta reference.");
    }
  } catch (err) {
    console.error(err);
    setFastaStatus("Couldn't load XIST.fasta automatically. Paste or upload a FASTA.", "#fca5a5");
    if (announce) {
      reportStatus("Failed to load XIST.fasta; paste or upload a FASTA.", true);
    }
  }
};

const handleCustomFasta = (label, contents) => {
  fastaInput.value = contents;
  setFastaStatus(`Loaded ${label}.`, "#22c55e");
  reportStatus(`Loaded ${label}.`);
};

const handleFileSelection = (file) => {
  if (!file) return;
  setFastaStatus(`Reading ${file.name}...`);
  const reader = new FileReader();
  reader.onload = () => {
    const text = typeof reader.result === "string" ? reader.result : "";
    if (!text.trim()) {
      setFastaStatus(`${file.name} appears empty.`, "#fca5a5");
      reportStatus(`${file.name} appears empty.`, true);
      return;
    }
    handleCustomFasta(file.name, text.trim());
  };
  reader.onerror = () => {
    setFastaStatus(`Unable to read ${file.name}.`, "#fca5a5");
    reportStatus(`Failed to read ${file.name}.`, true);
  };
  reader.readAsText(file);
};

let audioCtx;
let activeVoices = [];
let abortPlayback = null;

const midiToFreq = (midi) => 440 * Math.pow(2, (midi - 69) / 12);

const sanitizeSequence = (raw) => {
  if (!raw) return "";
  return raw
    .split(/\r?\n/)
    .filter((line) => !line.trim().startsWith(">"))
    .join("")
    .replace(/[^A-Za-z]/g, "")
    .toUpperCase();
};

const sequenceToNotes = (sequence, targetLength) => {
  if (!sequence.length) return [];
  const rootMidi = getRootMidi();
  const limit = Math.max(1, Math.min(targetLength || sequence.length, 2048));
  const letters = sequence.split("");
  const notes = [];
  for (let i = 0; i < limit; i += 1) {
    const letter = letters[i % letters.length];
    if (!letter) break;
    const mappedIndex = baseMap[letter];
    const scaleIndex = typeof mappedIndex === "number"
      ? mappedIndex
      : (letter.charCodeAt(0) + i) % MINOR_PENT_INTERVALS.length;
    const octaveLift = (Math.floor(i / MINOR_PENT_INTERVALS.length) % 4) * 12; // spread melody
    const midi = rootMidi + MINOR_PENT_INTERVALS[scaleIndex % MINOR_PENT_INTERVALS.length] + octaveLift;
    notes.push({ midi, base: letter, idx: i });
  }
  return notes;
};

const reportStatus = (message, isError = false) => {
  statusEl.textContent = message;
  statusEl.style.color = isError ? "#fca5a5" : "#cbd5f5";
};

const clampRelease = (lengthSec, releaseSec) => {
  if (releaseSec >= lengthSec) {
    return Math.max(0.01, lengthSec * 0.8);
  }
  return releaseSec;
};

const triggerNote = (ctx, midi, params) => {
  const {
    lengthSec,
    attackSec,
    decaySec,
    sustainLevel,
    releaseSec,
    filterCutoff,
    resonance,
    waveformType,
    lfoRate,
    lfoDepth,
    delayTime,
    delayFeedback,
    delayMix,
  } = params;

  const osc = ctx.createOscillator();
  const filter = ctx.createBiquadFilter();
  const gain = ctx.createGain();
  const delayNode = ctx.createDelay(2);
  const feedbackGain = ctx.createGain();
  const dryGain = ctx.createGain();
  const wetGain = ctx.createGain();
  let lfoOsc = null;
  let lfoGain = null;

  const start = ctx.currentTime;
  const limitedRelease = clampRelease(lengthSec, releaseSec);
  const sustainStart = Math.min(lengthSec, attackSec + decaySec);
  const releaseStart = Math.max(sustainStart, lengthSec - limitedRelease);

  osc.type = waveformType;
  osc.frequency.setValueAtTime(midiToFreq(midi), start);

  filter.type = "lowpass";
  filter.frequency.setValueAtTime(filterCutoff, start);
  filter.Q.setValueAtTime(resonance, start);

  const safeDelayMix = Math.min(1, Math.max(0, delayMix));
  dryGain.gain.setValueAtTime(1 - safeDelayMix, start);
  wetGain.gain.setValueAtTime(safeDelayMix, start);
  delayNode.delayTime.setValueAtTime(Math.min(2, Math.max(0, delayTime)), start);
  feedbackGain.gain.setValueAtTime(Math.min(0.95, Math.max(0, delayFeedback)), start);

  if (lfoRate > 0 && lfoDepth > 0) {
    lfoOsc = ctx.createOscillator();
    lfoGain = ctx.createGain();
    const depthHz = Math.max(0, filterCutoff * lfoDepth);
    lfoGain.gain.setValueAtTime(depthHz, start);
    lfoOsc.frequency.setValueAtTime(lfoRate, start);
    lfoOsc.connect(lfoGain);
    lfoGain.connect(filter.frequency);
    lfoOsc.start(start);
    lfoOsc.stop(start + lengthSec + limitedRelease + 0.5);
  }

  gain.gain.cancelScheduledValues(start);
  gain.gain.setValueAtTime(0.0001, start);
  gain.gain.linearRampToValueAtTime(1, start + attackSec);
  gain.gain.linearRampToValueAtTime(sustainLevel, start + sustainStart);
  gain.gain.setValueAtTime(sustainLevel, start + releaseStart);
  gain.gain.linearRampToValueAtTime(0.0001, start + releaseStart + limitedRelease);

  osc.connect(filter);
  filter.connect(gain);
  gain.connect(dryGain);
  dryGain.connect(ctx.destination);

  gain.connect(delayNode);
  delayNode.connect(wetGain);
  wetGain.connect(ctx.destination);
  delayNode.connect(feedbackGain);
  feedbackGain.connect(delayNode);

  const naturalStop = start + lengthSec + limitedRelease + 0.1;

  osc.start(start);
  osc.stop(naturalStop);

  const stop = () => {
    try {
      const now = ctx.currentTime;
      gain.gain.cancelScheduledValues(now);
      gain.gain.setValueAtTime(gain.gain.value, now);
      gain.gain.linearRampToValueAtTime(0.0001, now + 0.05);
      osc.stop(now + 0.05);
      if (lfoOsc) {
        lfoOsc.stop(now + 0.05);
        lfoOsc.disconnect();
      }
      if (lfoGain) {
        lfoGain.disconnect();
      }
      feedbackGain.gain.cancelScheduledValues(now);
      wetGain.gain.cancelScheduledValues(now);
      dryGain.gain.cancelScheduledValues(now);
      feedbackGain.gain.setValueAtTime(0, now);
      wetGain.gain.setValueAtTime(0, now);
      dryGain.gain.setValueAtTime(0, now);
      try { delayNode.disconnect(); } catch (e) { /* noop */ }
      try { feedbackGain.disconnect(); } catch (e) { /* noop */ }
      try { wetGain.disconnect(); } catch (e) { /* noop */ }
      try { dryGain.disconnect(); } catch (e) { /* noop */ }
      try { gain.disconnect(); } catch (e) { /* noop */ }
      try { filter.disconnect(); } catch (e) { /* noop */ }
      try { osc.disconnect(); } catch (e) { /* noop */ }
    } catch (err) {
      // already stopped
    }
  };

  return { stop };
};

const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

const stopActiveVoices = () => {
  activeVoices.forEach((voice) => voice.stop());
  activeVoices = [];
};

const stopPlaybackNow = (message = "Playback stopped.") => {
  if (abortPlayback) {
    abortPlayback.stopped = true;
    abortPlayback = null;
  }
  stopActiveVoices();
  reportStatus(message);
};

const playNotes = async (ctx, notes, abortHandle) => {
  for (let i = 0; i < notes.length; i += 1) {
    if (abortHandle.stopped) break;
    const params = gatherParams();
    params.releaseSec = clampRelease(params.lengthSec, params.releaseSec);
    const event = notes[i];
    const midi = clampToOctaveWindow(event.midi);
    pushHelixEvent(event.base, midi);
    const voice = triggerNote(ctx, midi, params);
    activeVoices.push(voice);
    await wait(params.lengthSec * 1000);
  }
};

const gatherParams = () => {
  const lengthMs = Number(controls.noteLength.value);
  return {
    lengthSec: Math.max(0.05, lengthMs / 1000),
    attackSec: Number(controls.attack.value),
    decaySec: Number(controls.decay.value),
    sustainLevel: Number(controls.sustain.value),
    releaseSec: Number(controls.release.value),
    filterCutoff: Number(controls.filterCutoff.value),
    resonance: Number(controls.resonance.value),
    waveformType: waveformSelect?.value || "sawtooth",
    lfoRate: Number(lfoRateInput?.value ?? 0),
    lfoDepth: Number(lfoDepthInput?.value ?? 0),
    delayTime: Math.max(0, Math.min(2, Number(controls.delayTime?.value ?? 0) / 1000)),
    delayFeedback: Math.max(0, Math.min(0.95, Number(controls.delayFeedback?.value ?? 0))),
    delayMix: Math.max(0, Math.min(1, Number(controls.delayMix?.value ?? 0))),
  };
};

const handlePlay = async () => {
  if (reminderActive) {
    reportStatus("Dismiss the silent-mode reminder to start playback.", true);
    return;
  }
  const cleaned = sanitizeSequence(fastaInput.value);
  if (!cleaned.length) {
    reportStatus("Please provide a FASTA sequence (letters A/C/G/T/U).", true);
    return;
  }

  const windowed = getSequenceWindow(cleaned);
  if (!windowed.sequence.length) {
    reportStatus("The selected sequence range contains no valid bases.", true);
    return;
  }

  const maxNotes = Number(maxNotesInput.value) || windowed.sequence.length;
  const notes = sequenceToNotes(windowed.sequence, maxNotes);
  if (!notes.length) {
    reportStatus("No valid bases found to convert into notes.", true);
    return;
  }

  if (!audioCtx) {
    audioCtx = new (window.AudioContext || window.webkitAudioContext)();
  }
  await audioCtx.resume();

  if (abortPlayback) {
    stopPlaybackNow("Restarting playback...");
  }

  resetHelix();

  abortPlayback = { stopped: false };
  const keyName = keySelect?.value || "A";
  const windowLabel = getOctaveRangeLabel();
  const baseRangeLabel = windowed.start && windowed.end
    ? `bases ${windowed.start}–${windowed.end}`
    : "entire sequence";
  reportStatus(`Playing ${notes.length} notes (${baseRangeLabel}) in ${keyName} minor pentatonic (${windowLabel}).`);

  try {
    await playNotes(audioCtx, notes, abortPlayback);
    if (!abortPlayback?.stopped) {
      reportStatus("Sequence complete. Try tweaking the controls for new colors.");
    }
  } catch (err) {
    console.error(err);
    reportStatus("Playback error. Check console for details.", true);
  } finally {
    abortPlayback = null;
    stopActiveVoices();
  }
};

playBtn.addEventListener("click", handlePlay);
stopBtn.addEventListener("click", () => stopPlaybackNow());

if (keySelect) {
  keySelect.addEventListener("change", () => {
    reportStatus(`Scale set to ${keySelect.value} minor pentatonic (${getOctaveRangeLabel()}).`);
  });
}

const handleOctaveChange = () => {
  reportStatus(`Octave window set to ${getOctaveRangeLabel()}.`);
};

[octaveMinSelect, octaveMaxSelect].forEach((select) => {
  select?.addEventListener("change", handleOctaveChange);
});

dismissReminderBtn?.addEventListener("click", dismissMobileReminder);

const handleSequenceWindowChange = () => {
  const { start, end } = normalizeSequenceWindowInputs();
  reportStatus(`Sequence window set to bases ${start}–${end}.`);
};

[sequenceStartInput, sequenceEndInput].forEach((input) => {
  input?.addEventListener("change", handleSequenceWindowChange);
});

reloadDefaultBtn?.addEventListener("click", () => loadDefaultFasta(true));

fastaFileInput?.addEventListener("change", (event) => {
  const file = event.target.files?.[0];
  if (file) {
    handleFileSelection(file);
  }
  event.target.value = "";
});

window.addEventListener("load", () => {
  initHelix();
  loadDefaultFasta(false);
  if (isMobileDevice()) {
    showMobileReminder();
  }
});
