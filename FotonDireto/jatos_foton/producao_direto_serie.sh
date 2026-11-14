#!/bin/bash
# Roda o SINAL nas 4 faixas, mas UMA DE CADA VEZ.
#
# O producao_direto.sh dispara as 4 em paralelo, que e' o jeito normal. Este
# aqui existe para quando a maquina ja' esta' ocupada com a producao do ruido:
# 4 processos do ROOT a mais nao cabem na memoria, e um estouro mataria a
# rodada longa. Em serie gasta ~450 MB e um nucleo so'.
#
#   ./producao_direto_serie.sh            100 mil eventos por faixa
#   ./producao_direto_serie.sh 5000       teste rapido

cd "$(dirname "$0")" || exit 1

N=${1:-100000}

echo "[$(date +%H:%M:%S)] SINAL em serie -- 4 faixas, $N eventos cada"

root -l -b -q "foton_direto.cc( 4.,  7., 43, $N, \"pthat04_07\")" > log_direto04_07.log 2>&1
echo "[$(date +%H:%M:%S)]   faixa 4-7 pronta"
root -l -b -q "foton_direto.cc( 7., 10., 44, $N, \"pthat07_10\")" > log_direto07_10.log 2>&1
echo "[$(date +%H:%M:%S)]   faixa 7-10 pronta"
root -l -b -q "foton_direto.cc(10., 15., 45, $N, \"pthat10_15\")" > log_direto10_15.log 2>&1
echo "[$(date +%H:%M:%S)]   faixa 10-15 pronta"
root -l -b -q "foton_direto.cc(15., -1., 46, $N, \"pthat15_inf\")" > log_direto15_inf.log 2>&1
echo "[$(date +%H:%M:%S)] sinal terminou"

echo
echo "=============== SECOES DE CHOQUE ==============="
grep -H "secao de choque" log_direto*.log | sed 's/log_//; s/.log:  secao de choque *:/ /'
