# Mojo and Rust Setup for Quality Per Position Benchmark

This document describes how to set up Mojo and Rust environments to run the Python vs Rust vs Mojo benchmark for the `quality_per_position_boxplot_data` algorithm.

**Algorithm spec:** All three implementations (Python, Rust, Mojo) must follow the same pseudo-code. See [quality-per-position-algorithm-spec.md](quality-per-position-algorithm-spec.md) for the canonical specification.

**Native binaries:** Rust and Mojo run as standalone executables (not Python extensions), so each can use language-specific optimizations. The benchmark harness invokes them via subprocess. Rust/Mojo read the FASTQ path from argv and write JSON to stdout.

## Mojo Setup (for Users New to Mojo)

### Platform support

Mojo provides wheels only for **Linux** (`manylinux_2_34_x86_64`, `manylinux_2_34_aarch64`) and **macOS** (`macosx_13_0_arm64`). There is **no Windows native support**. On Windows, use WSL2 or Docker.

---

### Windows: Use WSL2 (recommended)

Mojo does not run natively on Windows. Use WSL2 with Ubuntu 24.04 so you can run Mojo from Cursor and the terminal.

#### 1. Install WSL2 and Ubuntu

```powershell
wsl --install -d Ubuntu-24.04
```

Restart if prompted. Open Ubuntu from the Start menu.

#### 2. Open the project in WSL from Cursor

In Cursor: **File → Open Folder** (or `Ctrl+K Ctrl+O`), then choose the WSL path:

```
\\wsl$\Ubuntu-24.04\home\<user>\...\enpicom-quality-per-position-challenge
```

If the repo is under `C:\Users\...`, from WSL it is at `/mnt/c/Users/.../enpicom-quality-per-position-challenge`. Or clone the repo inside WSL (e.g. `~/Repos/...`) for better performance.

Or from a WSL terminal: `cd` to the project and run `cursor .` (if the Cursor CLI is installed in WSL).

This makes the integrated terminal run in Linux, so `mojo` will work.

#### 3. Create a Linux venv (do not reuse the Windows .venv)

The Windows `.venv` uses `win_amd64` wheels; Mojo has none. Create a new venv **inside WSL**:

```bash
# In WSL, from project root
uv venv
source .venv/bin/activate
uv pip install pytest   # for tests; omit if using uv sync from pyproject.toml
uv pip install mojo --extra-index-url https://modular.gateway.scarf.sh/simple/
mojo --version
```

You can keep both `.venv` (Windows) and a WSL-specific venv. Add `.venv` to `.gitignore` if it is not already; each platform will have its own.

#### 4. Optional: Docker alternative

If you prefer not to use WSL, run Mojo in a container. Run from WSL or a Linux shell (so `$PWD` works):

```bash
docker run -it --rm -v "$(pwd):/app" -w /app ubuntu:24.04 bash -c '
  apt-get update && apt-get install -y curl
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="$HOME/.local/bin:$PATH"
  uv venv && source .venv/bin/activate
  uv pip install mojo --extra-index-url https://modular.gateway.scarf.sh/simple/
  mojo --version
'
```

---

### Linux and macOS: Add Mojo to this project (uv)

1. Ensure the venv is active (or create one):

   ```bash
   uv venv && source .venv/bin/activate
   ```

2. Install Mojo with the Modular index:

   ```bash
   uv pip install mojo --extra-index-url https://modular.gateway.scarf.sh/simple/
   ```

3. Verify:

   ```bash
   mojo --version
   ```

### Linux and macOS: Add Mojo with Pixi (alternative)

If you prefer Pixi instead of uv:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
pixi init -c https://conda.modular.com/max-nightly/ -c conda-forge
pixi add mojo
pixi shell
mojo --version
```

**Note:** Use either uv or Pixi, not both.

### IDE Support

Install the **Mojo** extension from Modular for VS Code or Cursor:

- [VS Code Marketplace](https://marketplace.visualstudio.com/items?itemName=modular-mojotools.vscode-mojo)
- [Open VSX Registry](https://open-vsx.org/extension/modular-mojotools/vscode-mojo) (for Cursor and other editors)

### Quick Test

```bash
echo 'def main() raises: print("Hello, Mojo!")' > hello.mojo
mojo hello.mojo
# Expected: Hello, Mojo!
```

## Rust Setup

### Prerequisites

- [Rust toolchain](https://rustup.rs/) (rustc, cargo)
- Python 3.10+ (for maturin/PyO3 bindings)

### Installation

1. Install Rust:

   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

2. Build the Rust crate (when implemented):

   ```bash
   cd quality_per_position_rust
   maturin develop
   ```

## Running the Benchmark

1. Generate synthetic data (creates `tmp/` if needed):

   ```bash
   python scripts/generate_benchmark_fastq.py S
   ```

2. Run the benchmark harness:

   ```bash
   python scripts/benchmark_quality_per_position.py S
   ```

   For larger tiers (M, L, XL), replace `S` with the desired tier. Mojo and Rust columns show "—" if not installed.

## References

- [Mojo installation](https://docs.modular.com/mojo/manual/install)
- [Calling Mojo from Python](https://docs.modular.com/mojo/manual/python/mojo-from-python)
- [PyO3 (Rust–Python bindings)](https://pyo3.rs/)
