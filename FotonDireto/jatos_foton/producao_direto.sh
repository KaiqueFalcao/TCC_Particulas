#!/bin/bash
# Roda o SINAL (foton_direto.cc) nas MESMAS 4 faixas exclusivas de pThat
# usadas no ruido, para que a pureza possa ser calculada.
#
#   ./producao_direto.sh            100 mil eventos por faixa
#   ./producao_direto.sh 5000       teste rapido
#
# As sementes sao as mesmas do ruido (43, 44, 45, 46). Nao ha' problema:
# as duas amostras sao geradas por processos diferentes e nunca sao
# comparadas evento a evento, so' por secao de choque.

cd "$(dirname "$0")" || exit 1

N=${1:-100000}

echo "[$(date +%H:%M:%S)] SINAL -- 4 faixas exclusivas, $N eventos cada"
root -l -b -q "foton_direto.cc( 4.,  7., 43, $N, \"pthat04_07\")" > log_direto04_07.log 2>&1 &
root -l -b -q "foton_direto.cc( 7., 10., 44, $N, \"pthat07_10\")" > log_direto07_10.log 2>&1 &
root -l -b -q "foton_direto.cc(10., 15., 45, $N, \"pthat10_15\")" > log_direto10_15.log 2>&1 &
root -l -b -q "foton_direto.cc(15., -1., 46, $N, \"pthat15_inf\")" > log_direto15_inf.log 2>&1 &
wait
echo "[$(date +%H:%M:%S)] sinal terminou"

echo
echo "=============== SECOES DE CHOQUE ==============="
grep -H "secao de choque" log_direto*.log | sed 's/log_//; s/.log:  secao de choque *:/ /'
