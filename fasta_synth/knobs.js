const KNOB_ANGLE_MIN = -135;
const KNOB_ANGLE_MAX = 135;
const KNOB_ANGLE_RANGE = KNOB_ANGLE_MAX - KNOB_ANGLE_MIN;
const KNOB_DRAG_PIXELS_FOR_FULL_SWEEP = 180;

const setKnobAngle = (knobEl, pct) => {
  const dial = knobEl.querySelector(".knob__dial");
  if (!dial) return;
  const angle = KNOB_ANGLE_MIN + Math.min(1, Math.max(0, pct)) * KNOB_ANGLE_RANGE;
  dial.style.transform = `rotate(${angle}deg)`;
};

const bindRangeKnob = (knobEl, input) => {
  const min = Number(input.min);
  const max = Number(input.max);
  const step = Number(input.step) || (max - min) / 100 || 1;
  const defaultValue = Number(input.value);

  const render = () => {
    const val = Number(input.value);
    setKnobAngle(knobEl, (val - min) / (max - min));
  };

  let dragging = false;
  let startY = 0;
  let startVal = 0;

  const onPointerMove = (e) => {
    if (!dragging) return;
    const dy = startY - e.clientY;
    const delta = (dy / KNOB_DRAG_PIXELS_FOR_FULL_SWEEP) * (max - min);
    let next = startVal + delta;
    next = Math.round(next / step) * step;
    next = Math.min(max, Math.max(min, next));
    if (Number(input.value) !== next) {
      input.value = next;
      input.dispatchEvent(new Event("input", { bubbles: true }));
      render();
    }
  };

  const onPointerUp = () => {
    if (!dragging) return;
    dragging = false;
    knobEl.classList.remove("dragging");
    input.dispatchEvent(new Event("change", { bubbles: true }));
    document.removeEventListener("pointermove", onPointerMove);
    document.removeEventListener("pointerup", onPointerUp);
  };

  knobEl.addEventListener("pointerdown", (e) => {
    dragging = true;
    startY = e.clientY;
    startVal = Number(input.value);
    knobEl.classList.add("dragging");
    document.addEventListener("pointermove", onPointerMove);
    document.addEventListener("pointerup", onPointerUp);
  });

  knobEl.addEventListener("keydown", (e) => {
    const bump = (mult) => {
      const next = Math.min(max, Math.max(min, Number(input.value) + step * mult));
      input.value = next;
      input.dispatchEvent(new Event("input", { bubbles: true }));
      input.dispatchEvent(new Event("change", { bubbles: true }));
      render();
    };
    if (e.key === "ArrowUp" || e.key === "ArrowRight") { bump(1); e.preventDefault(); }
    if (e.key === "ArrowDown" || e.key === "ArrowLeft") { bump(-1); e.preventDefault(); }
  });

  knobEl.addEventListener("dblclick", () => {
    input.value = defaultValue;
    input.dispatchEvent(new Event("input", { bubbles: true }));
    input.dispatchEvent(new Event("change", { bubbles: true }));
    render();
  });

  render();
};

const bindSelectKnob = (knobEl, select) => {
  const valueOut = document.getElementById(`${select.id}Value`);
  const defaultIndex = select.selectedIndex;

  const render = () => {
    const count = select.options.length;
    const idx = select.selectedIndex;
    setKnobAngle(knobEl, count > 1 ? idx / (count - 1) : 0);
    if (valueOut) valueOut.textContent = select.options[idx]?.text ?? "";
  };

  let dragging = false;
  let startY = 0;
  let startIdx = 0;

  const onPointerMove = (e) => {
    if (!dragging) return;
    const dy = startY - e.clientY;
    const stepPx = 26;
    const steps = Math.round(dy / stepPx);
    const idx = Math.min(select.options.length - 1, Math.max(0, startIdx + steps));
    if (idx !== select.selectedIndex) {
      select.selectedIndex = idx;
      select.dispatchEvent(new Event("change", { bubbles: true }));
      render();
    }
  };

  const onPointerUp = () => {
    if (!dragging) return;
    dragging = false;
    knobEl.classList.remove("dragging");
    document.removeEventListener("pointermove", onPointerMove);
    document.removeEventListener("pointerup", onPointerUp);
  };

  knobEl.addEventListener("pointerdown", (e) => {
    dragging = true;
    startY = e.clientY;
    startIdx = select.selectedIndex;
    knobEl.classList.add("dragging");
    document.addEventListener("pointermove", onPointerMove);
    document.addEventListener("pointerup", onPointerUp);
  });

  knobEl.addEventListener("keydown", (e) => {
    const bump = (dir) => {
      const idx = Math.min(select.options.length - 1, Math.max(0, select.selectedIndex + dir));
      if (idx !== select.selectedIndex) {
        select.selectedIndex = idx;
        select.dispatchEvent(new Event("change", { bubbles: true }));
        render();
      }
    };
    if (e.key === "ArrowUp" || e.key === "ArrowRight") { bump(1); e.preventDefault(); }
    if (e.key === "ArrowDown" || e.key === "ArrowLeft") { bump(-1); e.preventDefault(); }
  });

  knobEl.addEventListener("dblclick", () => {
    select.selectedIndex = defaultIndex;
    select.dispatchEvent(new Event("change", { bubbles: true }));
    render();
  });

  render();
};

window.addEventListener("load", () => {
  document.querySelectorAll(".knob").forEach((knobEl) => {
    const targetId = knobEl.dataset.target;
    const type = knobEl.dataset.type;
    const target = document.getElementById(targetId);
    if (!target) return;
    if (type === "select") {
      bindSelectKnob(knobEl, target);
    } else {
      bindRangeKnob(knobEl, target);
    }
  });
});
