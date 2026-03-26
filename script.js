const MAX_SIZE = 6;
const EPS = 1e-12;

const sizeInput = document.getElementById("sizeInput");
const methodSelect = document.getElementById("methodSelect");
const generateBtn = document.getElementById("generateBtn");
const solveBtn = document.getElementById("solveBtn");
const matrixContainer = document.getElementById("matrixContainer");
const lOutput = document.getElementById("lOutput");
const uOutput = document.getElementById("uOutput");
const xOutput = document.getElementById("xOutput");
const errorBox = document.getElementById("errorBox");

function clampSize(value) {
  const parsed = Number(value);
  if (!Number.isInteger(parsed) || parsed < 1) return 1;
  if (parsed > MAX_SIZE) return MAX_SIZE;
  return parsed;
}

function clearOutputs() {
  lOutput.innerHTML = "";
  uOutput.innerHTML = "";
  xOutput.innerHTML = "";
}

function setError(message) {
  errorBox.textContent = message || "";
}

function generateInputs() {
  const n = clampSize(sizeInput.value);
  sizeInput.value = n;

  let html = '<table class="matrix-table"><thead><tr><th colspan="' + (n + 2) + '">Sistema Ax = b</th></tr></thead><tbody>';

  for (let i = 0; i < n; i += 1) {
    html += "<tr>";
    for (let j = 0; j < n; j += 1) {
      const defaultValue = i === j ? "1" : "0";
      html += '<td><input type="number" step="any" data-type="a" data-row="' + i + '" data-col="' + j + '" value="' + defaultValue + '"></td>';
    }
    html += '<td><span class="matrix-badge">|</span></td>';
    html += '<td><input type="number" step="any" data-type="b" data-row="' + i + '" value="0"></td>';
    html += "</tr>";
  }

  html += "</tbody></table>";
  matrixContainer.innerHTML = html;
  clearOutputs();
  setError("");
}

function readSystem() {
  const n = clampSize(sizeInput.value);
  const A = Array.from({ length: n }, () => Array(n).fill(0));
  const b = Array(n).fill(0);

  const aInputs = matrixContainer.querySelectorAll('input[data-type="a"]');
  const bInputs = matrixContainer.querySelectorAll('input[data-type="b"]');

  if (aInputs.length !== n * n || bInputs.length !== n) {
    throw new Error("Genera la matriz antes de resolver.");
  }

  for (const input of aInputs) {
    const i = Number(input.dataset.row);
    const j = Number(input.dataset.col);
    const value = Number(input.value);
    if (!Number.isFinite(value)) {
      throw new Error("Todos los valores de A deben ser numeros validos.");
    }
    A[i][j] = value;
  }

  for (const input of bInputs) {
    const i = Number(input.dataset.row);
    const value = Number(input.value);
    if (!Number.isFinite(value)) {
      throw new Error("Todos los valores de b deben ser numeros validos.");
    }
    b[i] = value;
  }

  return { A, b, n };
}

function doolittle(A) {
  const n = A.length;
  const L = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))
  );
  const U = Array.from({ length: n }, () => Array(n).fill(0));

  for (let k = 0; k < n; k += 1) {
    for (let j = k; j < n; j += 1) {
      let sum = 0;
      for (let s = 0; s < k; s += 1) {
        sum += L[k][s] * U[s][j];
      }
      U[k][j] = A[k][j] - sum;
    }

    if (Math.abs(U[k][k]) < EPS) {
      throw new Error("Doolittle fallo: pivote cero. Intenta otra matriz.");
    }

    for (let i = k + 1; i < n; i += 1) {
      let sum = 0;
      for (let s = 0; s < k; s += 1) {
        sum += L[i][s] * U[s][k];
      }
      L[i][k] = (A[i][k] - sum) / U[k][k];
    }
  }

  return { L, U };
}

function crout(A) {
  const n = A.length;
  const L = Array.from({ length: n }, () => Array(n).fill(0));
  const U = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))
  );

  for (let j = 0; j < n; j += 1) {
    for (let i = j; i < n; i += 1) {
      let sum = 0;
      for (let k = 0; k < j; k += 1) {
        sum += L[i][k] * U[k][j];
      }
      L[i][j] = A[i][j] - sum;
    }

    if (Math.abs(L[j][j]) < EPS) {
      throw new Error("Crout fallo: pivote cero. Intenta otra matriz.");
    }

    for (let i = j + 1; i < n; i += 1) {
      let sum = 0;
      for (let k = 0; k < j; k += 1) {
        sum += L[j][k] * U[k][i];
      }
      U[j][i] = (A[j][i] - sum) / L[j][j];
    }
  }

  return { L, U };
}

function cholesky(A) {
  const n = A.length;

  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
      if (Math.abs(A[i][j] - A[j][i]) > 1e-9) {
        throw new Error("Cholesky requiere una matriz simetrica.");
      }
    }
  }

  const L = Array.from({ length: n }, () => Array(n).fill(0));

  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j <= i; j += 1) {
      let sum = 0;
      for (let k = 0; k < j; k += 1) {
        sum += L[i][k] * L[j][k];
      }

      if (i === j) {
        const value = A[i][i] - sum;
        if (value <= EPS) {
          throw new Error("Cholesky requiere matriz definida positiva.");
        }
        L[i][j] = Math.sqrt(value);
      } else {
        if (Math.abs(L[j][j]) < EPS) {
          throw new Error("Cholesky fallo por division entre cero.");
        }
        L[i][j] = (A[i][j] - sum) / L[j][j];
      }
    }
  }

  const U = transpose(L);
  return { L, U };
}

function transpose(M) {
  return M[0].map((_, j) => M.map((row) => row[j]));
}

function forwardSubstitution(L, b) {
  const n = L.length;
  const y = Array(n).fill(0);

  for (let i = 0; i < n; i += 1) {
    let sum = 0;
    for (let j = 0; j < i; j += 1) {
      sum += L[i][j] * y[j];
    }

    if (Math.abs(L[i][i]) < EPS) {
      throw new Error("No se pudo resolver Ly = b (diagonal cero en L).");
    }

    y[i] = (b[i] - sum) / L[i][i];
  }

  return y;
}

function backwardSubstitution(U, y) {
  const n = U.length;
  const x = Array(n).fill(0);

  for (let i = n - 1; i >= 0; i -= 1) {
    let sum = 0;
    for (let j = i + 1; j < n; j += 1) {
      sum += U[i][j] * x[j];
    }

    if (Math.abs(U[i][i]) < EPS) {
      throw new Error("No se pudo resolver Ux = y (diagonal cero en U).");
    }

    x[i] = (y[i] - sum) / U[i][i];
  }

  return x;
}

function formatNumber(value) {
  if (Math.abs(value) < 1e-10) return "0";
  return Number(value.toFixed(6)).toString();
}

function renderMatrix(target, matrix, labelPrefix) {
  if (!matrix || matrix.length === 0) {
    target.innerHTML = "";
    return;
  }

  let html = '<table class="result-table"><tbody>';
  for (let i = 0; i < matrix.length; i += 1) {
    html += "<tr>";
    for (let j = 0; j < matrix[i].length; j += 1) {
      html += "<td>" + formatNumber(matrix[i][j]) + "</td>";
    }
    html += "</tr>";
  }
  html += "</tbody></table>";
  target.innerHTML = html;
}

function renderVector(target, vector, symbol) {
  let html = '<table class="result-table"><thead><tr><th>Variable</th><th>Valor</th></tr></thead><tbody>';
  for (let i = 0; i < vector.length; i += 1) {
    html += "<tr><td>" + symbol + (i + 1) + "</td><td>" + formatNumber(vector[i]) + "</td></tr>";
  }
  html += "</tbody></table>";
  target.innerHTML = html;
}

function solveSystem() {
  setError("");
  clearOutputs();

  try {
    const { A, b } = readSystem();
    const method = methodSelect.value;

    let decomp;
    if (method === "doolittle") {
      decomp = doolittle(A);
    } else if (method === "crout") {
      decomp = crout(A);
    } else {
      decomp = cholesky(A);
    }

    const y = forwardSubstitution(decomp.L, b);
    const x = backwardSubstitution(decomp.U, y);

    renderMatrix(lOutput, decomp.L, "l");
    renderMatrix(uOutput, decomp.U, "u");
    renderVector(xOutput, x, "x");
  } catch (error) {
    const message = error instanceof Error ? error.message : "Ocurrio un error inesperado.";
    setError(message);
  }
}

sizeInput.addEventListener("change", () => {
  sizeInput.value = clampSize(sizeInput.value);
});

generateBtn.addEventListener("click", generateInputs);
solveBtn.addEventListener("click", solveSystem);

generateInputs();
