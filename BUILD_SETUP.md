# Build Setup Guide for Windows

## Issue
The project currently cannot compile on Windows due to missing Visual Studio C++ Build Tools.

## Error Message
```
linking with `link.exe` failed: exit code: 1
note: in the Visual Studio installer, ensure the "C++ build tools" workload is selected
```

## Solution

### Option 1: Install Visual Studio Build Tools (Recommended)

1. Download [Visual Studio Build Tools](https://visualstudio.microsoft.com/downloads/)
2. Run the installer
3. Select **"Desktop development with C++"** workload
4. Complete the installation (requires ~7GB)
5. Restart your terminal/PowerShell
6. Run `cargo build` again

### Option 2: Install Full Visual Studio

1. Download [Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)
2. During installation, select **"Desktop development with C++"**
3. Complete installation
4. Restart terminal
5. Run `cargo build`

### Option 3: Use WSL (Windows Subsystem for Linux)

If you prefer a Linux environment:

```powershell
# In PowerShell (as Administrator)
wsl --install

# After restart, in WSL terminal:
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cd /mnt/c/Users/hemav/OneDrive/Desktop/bn254
cargo build
cargo test
```

## Verifying Installation

After installing the build tools, verify with:

```powershell
cargo --version
rustc --version
cargo build
```

## Next Steps

Once the build environment is set up:
1. `cargo build` - Compile the project
2. `cargo test` - Run all tests
3. `cargo test --lib` - Run only library tests
4. `cargo test --test integration` - Run integration tests
