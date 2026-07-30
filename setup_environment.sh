#!/bin/bash
set -e

echo "[+] Updating system packages and installing bioinformatics binaries..."
sudo apt-get update
sudo apt-get install -y default-jdk bcftools tabix wget unzip python3-pip

echo "[+] Creating local bin directory..."
mkdir -p bin

echo "[+] Downloading official PharmCAT v3.4.0..."
PHARMCAT_URL="https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat-3.4.0-all.jar"
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
elif [ -f "$HOME/pharmcat/vcf-preprocessor/requirements.txt" ]; then
    pip3 install --user -r "$HOME/pharmcat/vcf-preprocessor/requirements.txt"
fi

echo "[+] Automatically generating local SQLite genomic knowledgebase..."
if [ -f "init_db.py" ]; then
    python3 init_db.py
    echo "[✔] SQLite knowledgebase database generated successfully!"
else
    echo "[!] Warning: init_db.py not found in root directory. Skipping database generation."
fi

echo "[+] Configuring PATH and environment variables..."
WORKSPACE_DIR="$(pwd)"

# Avoid duplicate entries in .bashrc if run multiple times
grep -q "PHARMCAT_JAR" ~/.bashrc || echo "export PHARMCAT_JAR=\"$WORKSPACE_DIR/bin/pharmcat.jar\"" >> ~/.bashrc
grep -q "$WORKSPACE_DIR/bin" ~/.bashrc || echo "export PATH=\"\$PATH:$WORKSPACE_DIR/bin\"" >> ~/.bashrc

echo "[✔] Environment fully configured, database built, and ready!"
