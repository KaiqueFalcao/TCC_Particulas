#!/bin/bash
# Roda a producao do jatos_foton_v2.cc.
#
#   ./producao.sh                      producao padrao (ver abaixo)
#   ./producao.sh 5000 teto 5000       teste rapido, 5 mil em todas as faixas
#   ./producao.sh 1000000 teto         so' as faixas COM teto (as que se somam)
#   ./producao.sh 1000000 semteto      so' as faixas SEM teto (encaixadas)
#
# QUANTOS EVENTOS, E POR QUE DIFERENTE EM CADA FAIXA
# As tres primeiras faixas (4-7, 7-10, 10-15) levam N eventos cada. A faixa
# de cima (> 15) leva N_ALTO, dez vezes mais, porque e' ela que alimenta a
# regiao de pT alto do espectro -- justamente onde a pureza mais sobe e onde
# a estatistica era de poucos eventos por bin.
#
# Isso NAO desequilibra a soma: o peso de cada faixa e' sigma/N, com o N
# aceito lido do hist_sigma do proprio arquivo. Gerar dez vezes mais eventos
# na faixa de cima divide o peso dela por dez, e a secao de choque somada
# continua a mesma -- muda so' a barra de erro, que fica menor.
#
# Sao 4 rodadas por conjunto. A maquina tem poucos nucleos, entao os dois
# conjuntos vao em LEVAS SEPARADAS de 4, nunca 8 ao mesmo tempo.
#
# COM TETO  (4-7, 7-10, 10-15, >15): faixas exclusivas, podem ser SOMADAS.
# SEM TETO  (>4, >7, >10, >15): faixas encaixadas, servem para COMPARAR
#           entre si, mas somar contaria a regiao dura quatro vezes.

cd "$(dirname "$0")" || exit 1

N=${1:-1000000}
CONJ=${2:-teto}
N_ALTO=${3:-10000000}

roda_teto() {
  echo "[$(date +%H:%M:%S)] COM TETO -- 4 faixas exclusivas"
  echo "            faixas 4-7, 7-10 e 10-15 : $N eventos cada"
  echo "            faixa  > 15              : $N_ALTO eventos"
  root -l -b -q "jatos_foton_v2.cc(15., -1., 46, $N_ALTO, \"pthat15_inf\")" > log_pthat15_inf.log 2>&1 &
  root -l -b -q "jatos_foton_v2.cc( 4.,  7., 43, $N, \"pthat04_07\")" > log_pthat04_07.log 2>&1 &
  root -l -b -q "jatos_foton_v2.cc( 7., 10., 44, $N, \"pthat07_10\")" > log_pthat07_10.log 2>&1 &
  root -l -b -q "jatos_foton_v2.cc(10., 15., 45, $N, \"pthat10_15\")" > log_pthat10_15.log 2>&1 &
  wait
  echo "[$(date +%H:%M:%S)] com teto terminou"
}

roda_semteto() {
  echo "[$(date +%H:%M:%S)] SEM TETO -- 4 faixas encaixadas, $N eventos cada"
  root -l -b -q "jatos_foton_v2.cc( 4., -1., 43, $N)" > log_pthat04.log 2>&1 &
  root -l -b -q "jatos_foton_v2.cc( 7., -1., 44, $N)" > log_pthat07.log 2>&1 &
  root -l -b -q "jatos_foton_v2.cc(10., -1., 45, $N)" > log_pthat10.log 2>&1 &
  root -l -b -q "jatos_foton_v2.cc(15., -1., 46, $N)" > log_pthat15.log 2>&1 &
  wait
  echo "[$(date +%H:%M:%S)] sem teto terminou"
}

case "$CONJ" in
  teto)    roda_teto ;;
  semteto) roda_semteto ;;
  tudo)    roda_teto; roda_semteto ;;
  *) echo "conjunto desconhecido: $CONJ  (use tudo, teto ou semteto)"; exit 1 ;;
esac

echo
echo "=============== SECOES DE CHOQUE ==============="
grep -H "secao de choque" log_pthat*.log | sed 's/log_//; s/.log:  secao de choque *:/ /'
