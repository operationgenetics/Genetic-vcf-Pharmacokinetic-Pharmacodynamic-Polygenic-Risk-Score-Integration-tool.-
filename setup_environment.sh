#!/bin/bash
set -e

echo "[+] Updating system packages and installing bioinformatics binaries..."
sudo apt-get update
sudo apt-get install -y default-jdk bcftools tabix wget unzip

echo "[+] Creating local bin directory..."
mkdir -p bin

echo "[+] Downloading official PharmCAT..."
PHARMCAT_URL="https://github.com/pharmcat/pharmcat/releases/download/v3.0.0/pharmcat-3.0.0-all.jar"
if [ ! -f "bin/pharmcat.jar" ]; then
    wget -O bin/pharmcat.jar "$PHARMCAT_URL"
fi

echo "[+] Downloading PLINK2..."
if [ ! -f "bin/plink2" ]; then
    wget -O plink2_linux.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
    unzip -o plink2_linux.zip -d bin/
    chmod +x bin/plink2
    rm -f plink2_linux.zip LICENSE.txt
fi

echo "[+] Installing Python dependencies..."
if [ -f "requirements.txt" ]; then
    pip3 install --user -r requirements.txt
fi

echo "[+] Configuring PATH and environment variables..."
WORKSPACE_DIR="$(pwd)"
echo "export PHARMCAT_JAR=\"$WORKSPACE_DIR/bin/pharmcat.jar\"" >> ~/.bashrc
echo "export PATH=\"\$PATH:$WORKSPACE_DIR/bin\"" >> ~/.bashrc

echo "[✔] Environment fully configured and ready!"
echo "Run: source ~/.bashrc to apply PATH changes."
