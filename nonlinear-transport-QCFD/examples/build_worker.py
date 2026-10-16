"""Build the optional Windows worker from the bundled C# source."""
import os
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parent
compiler = Path(os.environ.get("WINDIR", "C:/Windows")) / "Microsoft.NET/Framework64/v4.0.30319/csc.exe"
if not compiler.is_file():
    raise SystemExit("The .NET Framework C# compiler was not found; use --pauli-action instead.")
subprocess.run([str(compiler), "/nologo", "/optimize+", "/platform:x64",
                "/out:" + str(root / "PauliProductWorker.exe"),
                str(root / "PauliProductWorker.cs")], check=True)
print("Compiled worker ready.")
