#!/bin/bash
set -e

echo "[+] Updating system packages and installing bioinformatics binaries..."
sudo apt-get update
sudo apt-get install -y openjdk-21-jdk-headless bcftools tabix

echo "[+] Downloading and installing PharmCAT..."
curl -fsSL https://get.pharmcat.org | bash

echo "[+] Installing PharmCAT Python dependencies..."
if [ -f "requirements.txt" ]; then
    pip3 install --user -r requirements.txt
elif [ -f "$HOME/pharmcat/vcf-preprocessor/requirements.txt" ]; then
    pip3 install --user -r "$HOME/pharmcat/vcf-preprocessor/requirements.txt"
fi

echo "[+] Configuring PATH..."
export PATH="$PATH:$HOME/pharmcat"
if ! grep -q "pharmcat" ~/.bashrc; then
    echo 'export PATH="$PATH:$HOME/pharmcat"' >> ~/.bashrc
fi

echo "[✔] Environment fully configured and ready!"
