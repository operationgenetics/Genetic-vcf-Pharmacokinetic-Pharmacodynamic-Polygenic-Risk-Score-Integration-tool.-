#!/bin/bash
set -euo pipefail

# Determine script directory
WORKSPACE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORKSPACE_DIR"

echo "[+] Updating system packages and installing bioinformatics binaries..."
SUDO=""
if [ "$(id -u)" -ne 0 ] && command -v sudo &> /dev/null; then
    SUDO="sudo"
fi

export DEBIAN_FRONTEND=noninteractive
$SUDO apt-get update -qq
$SUDO apt-get install -y default-jdk bcftools tabix wget unzip python3-pip python3-venv build-essential

echo "[+] Creating local bin directory..."
mkdir -p bin

echo "[+] Downloading official PharmCAT v3.4.0..."
PHARMCAT_URL="https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat-3.4.0-all.jar"
if [ ! -f "bin/pharmcat.jar" ]; then
    wget -q -O bin/pharmcat.jar "$PHARMCAT_URL"
fi

echo "[+] Downloading PLINK2..."
if [ ! -f "bin/plink2" ]; then
    wget -q -O plink2_linux.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
    unzip -q -o plink2_linux.zip plink2 -d bin/
    chmod +x bin/plink2
    rm -f plink2_linux.zip
fi

echo "[+] Setting up Python virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
fi

# Source virtualenv using absolute path
source "$WORKSPACE_DIR/venv/bin/activate"

echo "[+] Installing dependencies and packaging local module..."
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
pip install -e .

echo "[+] Automatically generating local SQLite genomic knowledgebase..."
python3 init_db.py

echo "[+] Configuring PATH and environment variables..."
TOUCH_BASHRC="${HOME}/.bashrc"
touch "$TOUCH_BASHRC"

grep -q "PHARMCAT_JAR" "$TOUCH_BASHRC" || echo "export PHARMCAT_JAR=\"$WORKSPACE_DIR/bin/pharmcat.jar\"" >> "$TOUCH_BASHRC"
grep -q "$WORKSPACE_DIR/bin" "$TOUCH_BASHRC" || echo "export PATH=\"\$PATH:$WORKSPACE_DIR/bin\"" >> "$TOUCH_BASHRC"

echo "[✔] Environment fully configured and database generated!"