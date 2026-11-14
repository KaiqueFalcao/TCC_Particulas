// ============================================================================
//  pureza.cc  --  quanto da amostra isolada e' sinal de verdade?
// ============================================================================
//
// Responde a pergunta central do trabalho:
//
//     de tudo que passa no corte de isolamento, que fracao e' foton DIRETO
//     (o que eu quero medir) e que fracao e' ruido (o que me engana)?
//
//     pureza(pT) = S(pT) / [ S(pT) + B(pT) ]
//
// com S e B em SECAO DE CHOQUE, nao em contagem. Tem que ser em secao de
// choque: sinal e ruido vem de amostras diferentes, geradas com numeros de
// eventos diferentes, e so' a secao de choque as poe na mesma escala fisica.
//
// Calcula a pureza em DOIS momentos:
//
//   ANTES  -- so' com o corte de isolamento
//   DEPOIS -- tambem com a marcacao por massa invariante, que joga fora o
//             foton que forma massa de pi0 ou de eta com algum parceiro
//             visivel
//
// ----------------------------------------------------------------------------
//  A MARCACAO AGORA E' APLICADA NAS DUAS AMOSTRAS
// ----------------------------------------------------------------------------
// Numa versao anterior deste codigo a marcacao so' era aplicada ao ruido, e o
// sinal entrava identico nas duas etapas da pureza. Aquilo era uma SUPOSICAO:
// que a marcacao nao remove foton direto nenhum.
//
// Agora a foton_direto.cc tambem marca, com os MESMOS criterios, e guarda o
// hist_pt_direto_pos. Entao a perda de sinal deixou de ser suposta e passou a
// ser MEDIDA. E' importante que seja: a pureza depois da marcacao e'
//
//     P = S_depois / (S_depois + B_depois)
//
// e um corte que remove 20% do fundo mas tambem 20% do sinal nao melhora
// nada. Com as duas amostras marcadas, a conta diz a verdade sozinha.
//
// O que a marcacao pega no sinal e' sempre ACASO: na amostra de sinal nao
// existe foton de decaimento, entao todo foton direto marcado formou massa de
// pi0 ou de eta com um foton qualquer do evento, sem ter nada a ver com ele.
//
// AS DUAS JANELAS: a marcacao procura parceiro na massa do pi0 (0,135) e na
// massa do eta (0,548). O eta decai em dois fotons como o pi0 e contamina do
// mesmo jeito; por isso ele tambem saiu da categoria "outro hadron" e foi
// juntado ao pi0 na categoria "pi0 + eta".
//
// ----------------------------------------------------------------------------
//  POR QUE AS DUAS AMOSTRAS PRECISAM DAS MESMAS FAIXAS DE pThat
// ----------------------------------------------------------------------------
// Uma rodada de pTHatMin = X e' INCOMPLETA para qualquer observavel perto de
// X: faltam os eventos de pThat abaixo de X que flutuaram para cima. Se o
// ruido vier de pThat > 4 e o sinal de pThat > 5, a razao entre os dois mistura
// essa incompletude com a fisica, e a pureza sai errada justamente na regiao
// de pT baixo, que e' onde ha' mais estatistica.
//
// Por isso este codigo exige as MESMAS 4 faixas exclusivas nas duas amostras,
// cada uma pesada pela sua propria secao de choque:
//
//     h_total = soma_i ( sigma_i / N_i ) * h_i
//
// ----------------------------------------------------------------------------
//  RODAR
// ----------------------------------------------------------------------------
//   root -l pureza.cc          // bins de 5 GeV/c na curva de pureza
//   root -l 'pureza.cc(1)'     // sem reagrupar, 1 GeV/c por bin (ruidoso)
//
// Precisa de 8 arquivos, as 4 faixas de cada amostra:
//   jatos_foton_pthat04_07.root ...   (ruido,  ./producao.sh)
//   foton_direto_pthat04_07.root ...  (sinal,  ./producao_direto.sh 100000)
//
// Saem 5 figuras:
//   pureza_espectros.png        sinal e ruido em secao de choque
//   pureza_espectros_norm.png   o mesmo normalizado por foton (area 1)
//   pureza_invariante.png       o mesmo em secao de choque invariante
//   pureza_eiso.png             E_iso das 4 origens, normalizado, com o corte
//   pureza_curva.png            a pureza em funcao do pT

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLine.h"

const Int_t NRODADAS = 4;
const Int_t NORIGENS = 3;   // 0 = pi0 + eta, 1 = outro hadron, 2 = fragmentacao

// Corte de isolamento, so' para desenhar a linha vertical no grafico de E_iso
const Double_t PZ_PT_ISO_MAX = 2.0;     // GeV/c

// ============================================================================
//  PARAMETROS DA SECAO DE CHOQUE INVARIANTE  --  conferir com o orientador
// ============================================================================
// A forma invariante do espectro e'
//
//     (1/2pi) d2sigma/(dpT dy)  =  sigma_bin / ( 2pi * largura_do_bin * dy )
//
// em mb/(GeV/c). O dy e' o intervalo de rapidez em que os fotons foram
// selecionados: o corte e' |eta| < 0,5 e, para foton (que nao tem massa),
// y = eta. Logo dy = 0,5 - (-0,5) = 1,0.
//
// Ligando PZ_USA_1_SOBRE_PT entra tambem o fator 1/pT, e a conta passa a ser
//
//     E d3sigma/dp3 = (1/(2pi pT)) d2sigma/(dpT dy)   em mb/(GeV/c)^2
//
// que e' a outra convencao da literatura. Tem que ser IGUAL ao que esta' no
// soma_pthat.cc, senao os dois graficos ficam em escalas diferentes.
const Double_t PZ_DY = 1.0;
const Bool_t   PZ_USA_1_SOBRE_PT = kFALSE;

const char *pz_sufixo[NRODADAS] = {"pthat04_07", "pthat07_10",
                                   "pthat10_15", "pthat15_inf"};

// Histogramas em escopo de arquivo: continuam acessiveis no prompt do ROOT
TH1F *hist_sinal     = 0;          // foton direto isolado, ANTES
TH1F *hist_sinal_pos = 0;          // o mesmo, DEPOIS da marcacao
TH1F *hist_fundo     = 0;          // soma das 3 origens, ANTES
TH1F *hist_fundo_pos = 0;          // o mesmo, DEPOIS
TH1F *hist_ruido    [NORIGENS];    // pi0+eta, outro hadron, frag, ANTES
TH1F *hist_ruido_pos[NORIGENS];    // o mesmo, DEPOIS

// E_iso: a soma de pT no cone, das 3 origens do ruido e do sinal
TH1F *hist_eiso    [NORIGENS];
TH1F *hist_eiso_sinal = 0;

TH1F *hist_pureza        = 0;      // S/(S+B) antes
TH1F *hist_pureza_pos    = 0;      // S/(S+B) depois
TH1F *hist_pureza_prompt = 0;      // antes, contando fragmentacao como sinal

TCanvas *p1 = 0;   // espectros em secao de choque
TCanvas *p2 = 0;   // curva de pureza
TCanvas *p3 = 0;   // E_iso das 4 origens
TCanvas *p4 = 0;   // espectros normalizados por foton
TCanvas *p5 = 0;   // espectros em secao de choque invariante

// Le um histograma das 4 faixas de uma amostra, cada uma multiplicada pelo
// peso sigma/N da sua faixa, e devolve a soma. O prefixo diz de que amostra:
// "jatos_foton" para o ruido, "foton_direto" para o sinal.
//
// Devolve 0 se faltar arquivo ou histograma, depois de dizer qual.
TH1F *soma_faixas(const char *prefixo, const char *nome_hist,
                  const char *apelido)
{
  TH1F *total = 0;
  char arq[200];

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(arq, "%s_%s.root", prefixo, pz_sufixo[i]);

    TFile *f = TFile::Open(arq);
    if (!f || f->IsZombie())
    {
      printf("  ERRO: nao consegui abrir %s\n", arq);
      if (total) delete total;
      return 0;
    }

    TH1D *info = (TH1D*) f->Get("hist_sigma");
    TH1F *h    = (TH1F*) f->Get(nome_hist);

    if (!info || !h)
    {
      printf("  ERRO: %s nao tem %s\n", arq, info ? nome_hist : "hist_sigma");
      f->Close();
      if (total) delete total;
      return 0;
    }

    Double_t peso = info->GetBinContent(5);   // bin 5 = sigma / n_aceitos

    // O clone nasce com nome temporario e so' depois e' rebatizado: clonar
    // ja' com o nome final o criaria dentro do TFile, com nome repetido.
    TH1F *copia = (TH1F*) h->Clone("pz_tmp");
    copia->SetDirectory(0);          // solta do arquivo, sobrevive ao Close
    copia->Scale(peso);              // contagem -> mb

    if (!total) { total = copia; total->SetName(apelido); }
    else        { total->Add(copia); delete copia; }

    f->Close();
  }
  return total;
}

// Divide um espectro pela sua propria area, de modo que a area debaixo da
// curva passe a valer 1. O eixo deixa de ser secao de choque e passa a ser
// (1/N) dN/dX -- so' a FORMA da distribuicao.
//
// Scale(c, "width") divide cada bin pela largura e multiplica por c; com
// c = 1/integral sai a densidade de area 1.
//
// ATENCAO: curva normalizada NAO serve para calcular pureza. A pureza e'
// razao de secoes de choque, e a normalizacao joga fora justamente a
// informacao de tamanho que distingue sinal de fundo.
TH1F *faz_normalizado(TH1F *h, const char *nome)
{
  TH1F *nm = (TH1F*) h->Clone(nome);
  nm->SetDirectory(0);               // copia de grafico: fora do .root

  Double_t area = nm->Integral();    // soma dos bins visiveis
  if (area > 0) nm->Scale(1. / area, "width");

  nm->SetStats(0);
  return nm;
}

// Converte um espectro de "sigma [mb] por bin" para SECAO DE CHOQUE
// INVARIANTE, dividindo cada bin por 2pi * largura * dy (e por pT, se a
// chave PZ_USA_1_SOBRE_PT estiver ligada).
TH1F *faz_invariante(TH1F *h, const char *nome)
{
  TH1F *inv = (TH1F*) h->Clone(nome);
  inv->SetDirectory(0);              // copia de grafico: fora do .root

  for (Int_t b = 1; b <= inv->GetNbinsX(); b++)
  {
    Double_t fator = 2. * TMath::Pi() * inv->GetBinWidth(b) * PZ_DY;
    if (PZ_USA_1_SOBRE_PT) fator *= inv->GetBinCenter(b);
    if (fator <= 0) continue;

    inv->SetBinContent(b, inv->GetBinContent(b) / fator);
    inv->SetBinError  (b, inv->GetBinError(b)   / fator);
  }

  inv->SetYTitle(PZ_USA_1_SOBRE_PT
      ? "E d^{3}#sigma/dp^{3}  [mb/(GeV/c)^{2}]"
      : "(1/2#pi) d^{2}#sigma/(dp_{T}dy)  [mb/(GeV/c)]");
  inv->SetStats(0);
  return inv;
}

// Preenche "saida" com a pureza bin a bin, a partir das quatro componentes.
//
//     P = X / (X + Y),  X = sinal,  Y = fundo
//
// Os erros sao propagados na mao:
//
//     dP/dX =  Y / (X+Y)^2        dP/dY = -X / (X+Y)^2
//     sigma_P^2 = [ Y^2 sigma_X^2 + X^2 sigma_Y^2 ] / (X+Y)^4
//
// NAO da' para usar TH1::Divide com a opcao "B" (binomial): binomial vale
// quando o numerador e' um subconjunto do denominador medido na MESMA
// amostra. Aqui sinal e fundo vem de amostras independentes.
//
// frag_e_sinal muda de convencao: se for verdadeiro, a fragmentacao conta
// como sinal em vez de fundo.
void calcula_pureza(TH1F *S, TH1F *Bp, TH1F *Bo, TH1F *Bf,
                    TH1F *saida, Bool_t frag_e_sinal)
{
  for (Int_t b = 1; b <= S->GetNbinsX(); b++)
  {
    Double_t s  = S ->GetBinContent(b), es  = S ->GetBinError(b);
    Double_t bp = Bp->GetBinContent(b), ebp = Bp->GetBinError(b);
    Double_t bo = Bo->GetBinContent(b), ebo = Bo->GetBinError(b);
    Double_t bf = Bf->GetBinContent(b), ebf = Bf->GetBinError(b);

    Double_t tot = s + bp + bo + bf;
    if (tot <= 0) continue;              // bin vazio nas duas amostras

    Double_t X, Y, eX2, eY2;

    if (frag_e_sinal)
    {
      X = s + bf;   eX2 = es*es + ebf*ebf;
      Y = bp + bo;  eY2 = ebp*ebp + ebo*ebo;
    }
    else
    {
      X = s;             eX2 = es*es;
      Y = bp + bo + bf;  eY2 = ebp*ebp + ebo*ebo + ebf*ebf;
    }

    // SE O FUNDO NAO TEM NENHUM EVENTO DE MONTE CARLO NESTE BIN, a pureza
    // nao e' 100%: e' DESCONHECIDA. Zero eventos gerados nao quer dizer
    // zero secao de choque, quer dizer que a amostra e' pequena demais para
    // dizer alguma coisa aqui. O ponto fica de fora do grafico.
    if (Y <= 0) continue;

    saida->SetBinContent(b, X / tot);
    saida->SetBinError  (b, TMath::Sqrt(Y*Y*eX2 + X*X*eY2) / (tot*tot));
  }
}

// rebin : quantos bins de 1 GeV/c agrupar na curva de pureza. O espectro
//         e' desenhado sempre na granularidade original; so' a pureza e'
//         reagrupada, porque ela e' uma razao e fica ilegivel com pouca
//         estatistica por bin.
void pureza(Int_t rebin = 5)
{
  const char *curto[NORIGENS] = {"pi0eta", "outro", "frag"};
  const char *nome_org[NORIGENS] = {"#pi^{0} + #eta", "outro hadron",
                                    "fragmentacao"};
  char nome_hist[100], apelido[100];

  printf("\n================= PUREZA =================\n");

  // ================================================================
  //  1) SINAL
  // ================================================================
  // Dois histogramas: antes e depois da marcacao por massa invariante. O
  // sinal passa pelo MESMO corte que o ruido, entao a perda de sinal e'
  // medida, nao suposta.
  hist_sinal     = soma_faixas("foton_direto", "hist_pt_direto_iso",
                               "hist_sinal");
  hist_sinal_pos = soma_faixas("foton_direto", "hist_pt_direto_pos",
                               "hist_sinal_pos");
  hist_eiso_sinal = soma_faixas("foton_direto", "hist_eiso_direto",
                                "hist_eiso_sinal");

  if (!hist_sinal || !hist_sinal_pos || !hist_eiso_sinal)
  {
    printf("\n  Falta a producao do SINAL nas 4 faixas de pThat, ou ela e'\n");
    printf("  anterior a marcacao por massa no sinal. Rode:\n\n");
    printf("    ./producao_direto.sh 100000\n\n");
    return;
  }

  // ================================================================
  //  2) RUIDO, antes e depois
  // ================================================================
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome_hist, "hist_pt_iso_%s", curto[o]);
    sprintf(apelido,   "hist_ruido_%s",  curto[o]);
    hist_ruido[o] = soma_faixas("jatos_foton", nome_hist, apelido);

    sprintf(nome_hist, "hist_pt_pos_%s",    curto[o]);
    sprintf(apelido,   "hist_ruido_pos_%s", curto[o]);
    hist_ruido_pos[o] = soma_faixas("jatos_foton", nome_hist, apelido);

    sprintf(nome_hist, "hist_eiso_%s", curto[o]);
    sprintf(apelido,   "hist_eiso_r_%s", curto[o]);
    hist_eiso[o] = soma_faixas("jatos_foton", nome_hist, apelido);

    if (!hist_ruido[o] || !hist_ruido_pos[o] || !hist_eiso[o])
    {
      printf("\n  Falta a producao do RUIDO, ou ela e' anterior a juncao\n");
      printf("  pi0 + eta, que mudou o nome dos histogramas.\n");
      printf("  Rode:\n    ./producao.sh\n");
      return;
    }
  }

  hist_fundo = (TH1F*) hist_ruido[0]->Clone("hist_fundo");
  hist_fundo->SetDirectory(0);
  hist_fundo->Add(hist_ruido[1]);
  hist_fundo->Add(hist_ruido[2]);

  hist_fundo_pos = (TH1F*) hist_ruido_pos[0]->Clone("hist_fundo_pos");
  hist_fundo_pos->SetDirectory(0);
  hist_fundo_pos->Add(hist_ruido_pos[1]);
  hist_fundo_pos->Add(hist_ruido_pos[2]);

  // ================================================================
  //  3) PUREZA
  // ================================================================
  // Copias reagrupadas, usadas so' na curva. O espectro do painel 1
  // continua na granularidade original, de 1 GeV/c.
  TH1F *S  = (TH1F*) hist_sinal    ->Clone("pz_S");   S ->SetDirectory(0);
  TH1F *Bp = (TH1F*) hist_ruido[0] ->Clone("pz_Bp");  Bp->SetDirectory(0);
  TH1F *Bo = (TH1F*) hist_ruido[1] ->Clone("pz_Bo");  Bo->SetDirectory(0);
  TH1F *Bf = (TH1F*) hist_ruido[2] ->Clone("pz_Bf");  Bf->SetDirectory(0);

  // Na etapa "depois", o sinal tambem passou pela marcacao.
  TH1F *Sd  = (TH1F*) hist_sinal_pos   ->Clone("pz_Sd");  Sd ->SetDirectory(0);
  TH1F *Bpd = (TH1F*) hist_ruido_pos[0]->Clone("pz_Bpd"); Bpd->SetDirectory(0);
  TH1F *Bod = (TH1F*) hist_ruido_pos[1]->Clone("pz_Bod"); Bod->SetDirectory(0);
  TH1F *Bfd = (TH1F*) hist_ruido_pos[2]->Clone("pz_Bfd"); Bfd->SetDirectory(0);

  if (rebin > 1)
  {
    S->Rebin(rebin);   Bp ->Rebin(rebin); Bo ->Rebin(rebin); Bf ->Rebin(rebin);
    Sd->Rebin(rebin);  Bpd->Rebin(rebin); Bod->Rebin(rebin); Bfd->Rebin(rebin);
  }

  hist_pureza = (TH1F*) S->Clone("hist_pureza");
  hist_pureza->SetDirectory(0);  hist_pureza->Reset();

  hist_pureza_pos = (TH1F*) S->Clone("hist_pureza_pos");
  hist_pureza_pos->SetDirectory(0);  hist_pureza_pos->Reset();

  hist_pureza_prompt = (TH1F*) S->Clone("hist_pureza_prompt");
  hist_pureza_prompt->SetDirectory(0);  hist_pureza_prompt->Reset();

  calcula_pureza(S,  Bp,  Bo,  Bf,  hist_pureza,        kFALSE);
  calcula_pureza(Sd, Bpd, Bod, Bfd, hist_pureza_pos,    kFALSE);
  calcula_pureza(S,  Bp,  Bo,  Bf,  hist_pureza_prompt, kTRUE);

  // ================================================================
  //  4) GRAFICOS
  // ================================================================
  // --- p1: os espectros, sinal e ruido no mesmo eixo
  p1 = new TCanvas("p1", "Sinal e ruido", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);

  Int_t cor_r[NORIGENS] = {kBlue, kGreen + 2, kRed};

  hist_fundo->SetTitle("#gamma isolados: sinal e ruido");
  hist_fundo->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  hist_fundo->SetYTitle("#sigma [mb] por bin");
  hist_fundo->SetStats(0);
  hist_fundo->SetLineColor(kBlack);
  hist_fundo->SetLineWidth(3);

  // O piso do eixo tem que caber o MENOR bin nao vazio das duas amostras.
  // Sem isso o ROOT escolhe o piso pelo histograma desenhado primeiro (o
  // fundo) e a cauda do sinal, que vai mais fundo, some embaixo do eixo.
  {
    Double_t piso = 1e30;
    for (Int_t b = 1; b <= hist_fundo->GetNbinsX(); b++)
    {
      Double_t cs = hist_sinal->GetBinContent(b);
      Double_t cb = hist_fundo->GetBinContent(b);
      if (cs > 0 && cs < piso) piso = cs;
      if (cb > 0 && cb < piso) piso = cb;
    }
    if (piso < 1e29) hist_fundo->SetMinimum(piso * 0.4);
  }

  hist_fundo->Draw("HIST E");

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    hist_ruido[o]->SetLineColor(cor_r[o]);
    hist_ruido[o]->SetLineWidth(1);
    hist_ruido[o]->SetLineStyle(2);
    hist_ruido[o]->Draw("HIST E SAME");
  }

  hist_sinal->SetLineColor(kMagenta + 2);
  hist_sinal->SetLineWidth(3);
  hist_sinal->Draw("HIST E SAME");

  TLegend *l1 = new TLegend(0.50, 0.63, 0.89, 0.89);
  l1->SetTextSize(0.030);
  l1->AddEntry(hist_sinal,    "SINAL: #gamma direto", "l");
  l1->AddEntry(hist_fundo,    "RUIDO: soma das 3",    "l");
  l1->AddEntry(hist_ruido[0], "ruido: #pi^{0} + #eta", "l");
  l1->AddEntry(hist_ruido[1], "ruido: outro hadron",   "l");
  l1->AddEntry(hist_ruido[2], "ruido: fragmentacao",   "l");
  l1->Draw();

  p1->Update();
  p1->SaveAs("pureza_espectros.png");

  // --- p2: a curva de pureza, antes e depois da marcacao
  p2 = new TCanvas("p2", "Pureza", 900, 700);
  gPad->SetLeftMargin(0.16);
  gPad->SetGridy();

  hist_pureza->SetTitle("pureza da amostra de #gamma isolados");
  hist_pureza->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  hist_pureza->SetYTitle("pureza = S / (S + B)");
  hist_pureza->SetStats(0);
  hist_pureza->SetMinimum(0.);
  hist_pureza->SetMaximum(1.05);
  hist_pureza->SetLineColor(kMagenta + 2);
  hist_pureza->SetLineWidth(3);
  hist_pureza->SetMarkerStyle(20);
  hist_pureza->SetMarkerColor(kMagenta + 2);
  hist_pureza->Draw("E1");

  hist_pureza_pos->SetLineColor(kAzure + 2);
  hist_pureza_pos->SetLineWidth(3);
  hist_pureza_pos->SetMarkerStyle(21);
  hist_pureza_pos->SetMarkerColor(kAzure + 2);
  hist_pureza_pos->Draw("E1 SAME");

  TLegend *l2 = new TLegend(0.40, 0.76, 0.89, 0.89);
  l2->SetTextSize(0.030);
  l2->AddEntry(hist_pureza,     "so' isolamento",          "lp");
  l2->AddEntry(hist_pureza_pos, "isolamento + massa m_{#gamma#gamma}", "lp");
  l2->SetHeader("marcacao aplicada nas duas amostras");
  l2->Draw();

  p2->Update();
  p2->SaveAs("pureza_curva.png");

  // --- p3: E_iso das QUATRO origens, cada curva com area 1.
  // E' o grafico que explica o isolamento: mostra, para cada tipo de foton,
  // como a atividade em volta dele se distribui. O foton direto nasce
  // sozinho e empilha perto de zero; o de decaimento vem dentro de um jato
  // e espalha para a direita. A linha tracejada e' o corte de 2 GeV/c: a
  // fracao de cada curva a esquerda dela e' a eficiencia do isolamento
  // para aquele tipo.
  p3 = new TCanvas("p3", "E_iso por origem", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  {
    // sinal primeiro, depois as tres origens do ruido
    TH1F *bruto_e[4] = {hist_eiso_sinal, hist_eiso[0],
                        hist_eiso[1],    hist_eiso[2]};
    Int_t cor_e[4] = {kMagenta + 2, kBlue, kGreen + 2, kRed};
    const char *nome_e[4] = {"SINAL: #gamma direto", "#pi^{0} + #eta",
                             "outro hadron", "fragmentacao"};
    const char *ap_e[4] = {"pz_eiso_n_sinal", "pz_eiso_n_pi0eta",
                           "pz_eiso_n_outro", "pz_eiso_n_frag"};

    TH1F *norm_e[4];
    Double_t topo_e = 0, piso_e = 1e30;

    for (Int_t i = 0; i < 4; i++)
    {
      norm_e[i] = faz_normalizado(bruto_e[i], ap_e[i]);
      if (norm_e[i]->GetMaximum() > topo_e) topo_e = norm_e[i]->GetMaximum();
      for (Int_t b = 1; b <= norm_e[i]->GetNbinsX(); b++)
      {
        Double_t v = norm_e[i]->GetBinContent(b);
        if (v > 0 && v < piso_e) piso_e = v;
      }
    }
    if (piso_e > 1e29) piso_e = 1e-6;

    TLegend *l3 = new TLegend(0.46, 0.66, 0.89, 0.89);
    l3->SetTextSize(0.031);

    for (Int_t i = 0; i < 4; i++)
    {
      norm_e[i]->SetTitle("isolamento: sinal e as tres origens do ruido");
      norm_e[i]->SetXTitle("E_{iso} = #Sigma pT^{cone} (R < 0,4) [GeV/c]");
      norm_e[i]->SetYTitle("(1/N_{#gamma}) dN/dE_{iso}  [c/GeV]");
      norm_e[i]->SetLineColor(cor_e[i]);
      norm_e[i]->SetLineWidth(i == 0 ? 3 : 2);
      norm_e[i]->SetMinimum(piso_e * 0.5);
      norm_e[i]->SetMaximum(topo_e * 3.);
      norm_e[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
      l3->AddEntry(norm_e[i], nome_e[i], "l");
    }

    TLine *lin3 = new TLine(PZ_PT_ISO_MAX, piso_e * 0.5,
                            PZ_PT_ISO_MAX, topo_e * 3.);
    lin3->SetLineColor(kBlack);
    lin3->SetLineStyle(2);
    lin3->SetLineWidth(2);
    lin3->Draw();
    l3->AddEntry(lin3, "corte de isolamento", "l");
    l3->Draw();

    // ---- tabela: eficiencia do isolamento das quatro
    printf("  EFICIENCIA DO ISOLAMENTO (fracao com E_iso < %.1f GeV/c)\n",
           PZ_PT_ISO_MAX);
    Int_t b_corte = bruto_e[0]->FindBin(PZ_PT_ISO_MAX) - 1;
    for (Int_t i = 0; i < 4; i++)
    {
      Int_t nbe = bruto_e[i]->GetNbinsX();
      Double_t tudo   = bruto_e[i]->Integral(0, nbe + 1);
      Double_t aceito = bruto_e[i]->Integral(1, b_corte);
      printf("    %-22s %8.2f%%   (<E_iso> = %6.3f GeV/c)\n",
             i == 0 ? "sinal (direto)" : curto[i - 1],
             tudo ? 100. * aceito / tudo : 0., bruto_e[i]->GetMean());
    }
    printf("------------------------------------------\n");
  }
  p3->Update();
  p3->SaveAs("pureza_eiso.png");

  // --- p4: os mesmos espectros de pT, normalizados por foton (area 1).
  // Aqui NAO se le pureza nenhuma -- as alturas relativas foram jogadas
  // fora de proposito. O que se le e' a FORMA: o sinal cai mais devagar
  // em pT do que o fundo de decaimento, e e' por isso que a pureza sobe
  // com o pT.
  p4 = new TCanvas("p4", "Espectros normalizados", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  {
    TH1F *bruto_n[4] = {hist_sinal, hist_ruido[0], hist_ruido[1],
                        hist_ruido[2]};
    Int_t cor_n[4] = {kMagenta + 2, kBlue, kGreen + 2, kRed};
    const char *nome_n[4] = {"SINAL: #gamma direto", "ruido: #pi^{0} + #eta",
                             "ruido: outro hadron", "ruido: fragmentacao"};
    const char *ap_n[4] = {"pz_pt_n_sinal", "pz_pt_n_pi0eta",
                           "pz_pt_n_outro", "pz_pt_n_frag"};

    TLegend *l4 = new TLegend(0.46, 0.70, 0.89, 0.89);
    l4->SetTextSize(0.031);

    for (Int_t i = 0; i < 4; i++)
    {
      TH1F *nm = faz_normalizado(bruto_n[i], ap_n[i]);
      nm->SetTitle("#gamma isolados, normalizado por #gamma");
      nm->SetXTitle("p_{T}^{#gamma} [GeV/c]");
      nm->SetYTitle("(1/N_{#gamma}) dN/dp_{T}  [c/GeV]");
      nm->SetLineColor(cor_n[i]);
      nm->SetLineWidth(i == 0 ? 3 : 2);
      nm->Draw(i == 0 ? "HIST E" : "HIST E SAME");
      l4->AddEntry(nm, nome_n[i], "l");
    }
    l4->SetHeader("area 1 em cada curva -- nao e' pureza");
    l4->Draw();
  }
  p4->Update();
  p4->SaveAs("pureza_espectros_norm.png");

  // --- p5: SECAO DE CHOQUE INVARIANTE do sinal e do ruido, antes e depois
  // da marcacao. Mesma fisica do p1, outra escala: cada bin dividido por
  // 2pi * largura * dy. E' a forma que se compara com dado publicado,
  // porque nao depende de quao largo o bin foi escolhido.
  p5 = new TCanvas("p5", "Secao de choque invariante", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.17);
  {
    TH1F *inv_s  = faz_invariante(hist_sinal,     "pz_inv_sinal");
    TH1F *inv_sd = faz_invariante(hist_sinal_pos, "pz_inv_sinal_pos");
    TH1F *inv_b  = faz_invariante(hist_fundo,     "pz_inv_fundo");
    TH1F *inv_bd = faz_invariante(hist_fundo_pos, "pz_inv_fundo_pos");

    // o piso do eixo tem que caber o menor bin nao vazio das quatro curvas
    Double_t piso = 1e30;
    TH1F *todas[4] = {inv_s, inv_sd, inv_b, inv_bd};
    for (Int_t i = 0; i < 4; i++)
      for (Int_t b = 1; b <= todas[i]->GetNbinsX(); b++)
      {
        Double_t v = todas[i]->GetBinContent(b);
        if (v > 0 && v < piso) piso = v;
      }
    if (piso > 1e29) piso = 1e-12;

    inv_b->SetTitle("#gamma isolados -- secao de choque invariante");
    inv_b->SetXTitle("p_{T}^{#gamma} [GeV/c]");
    inv_b->SetLineColor(kBlack);
    inv_b->SetLineWidth(3);
    inv_b->SetMinimum(piso * 0.4);
    inv_b->Draw("HIST E");

    inv_bd->SetLineColor(kBlack);
    inv_bd->SetLineWidth(2);
    inv_bd->SetLineStyle(2);
    inv_bd->Draw("HIST E SAME");

    inv_s->SetLineColor(kMagenta + 2);
    inv_s->SetLineWidth(3);
    inv_s->Draw("HIST E SAME");

    inv_sd->SetLineColor(kMagenta + 2);
    inv_sd->SetLineWidth(2);
    inv_sd->SetLineStyle(2);
    inv_sd->Draw("HIST E SAME");

    TLegend *l5 = new TLegend(0.40, 0.68, 0.89, 0.89);
    l5->SetTextSize(0.030);
    l5->AddEntry(inv_s,  "sinal, so' isolamento",  "l");
    l5->AddEntry(inv_sd, "sinal, + massa",         "l");
    l5->AddEntry(inv_b,  "ruido, so' isolamento",  "l");
    l5->AddEntry(inv_bd, "ruido, + massa",         "l");

    char cab[80];
    sprintf(cab, "dy = %.1f%s", PZ_DY,
            PZ_USA_1_SOBRE_PT ? ", com 1/p_{T}" : ", sem 1/p_{T}");
    l5->SetHeader(cab);
    l5->Draw();
  }
  p5->Update();
  p5->SaveAs("pureza_invariante.png");

  // ================================================================
  //  5) TABELA
  // ================================================================
  printf("  Sinal  : foton direto isolado, PromptPhoton:all\n");
  printf("  Ruido  : (pi0 + eta) + outro hadron + fragmentacao, HardQCD:all\n");
  printf("  As duas amostras somam as MESMAS 4 faixas de pThat, cada uma\n");
  printf("  pesada pela sua propria secao de choque, e passam pelo MESMO\n");
  printf("  corte de massa invariante, nas duas janelas (pi0 e eta).\n");
  printf("------------------------------------------\n");
  printf("  %10s %11s %11s %7s %19s %19s\n", "pT [GeV/c]",
         "S [mb]", "B [mb]", "N_B", "pureza so' isol.", "pureza + massa");
  printf("  (N_B = numero EFETIVO de eventos MC de fundo no bin, (c/e)^2)\n");

  for (Int_t b = 1; b <= S->GetNbinsX(); b++)
  {
    Double_t s  = S->GetBinContent(b);
    Double_t bg = Bp->GetBinContent(b) + Bo->GetBinContent(b)
                + Bf->GetBinContent(b);
    if (s + bg <= 0) continue;

    // Numero efetivo de eventos num histograma com pesos:
    //    n_ef = (soma w)^2 / (soma w^2) = (conteudo / erro)^2
    // E' quantos eventos nao pesados dariam a mesma incerteza relativa.
    Double_t eb = TMath::Sqrt(Bp->GetBinError(b)*Bp->GetBinError(b)
                            + Bo->GetBinError(b)*Bo->GetBinError(b)
                            + Bf->GetBinError(b)*Bf->GetBinError(b));
    Double_t nB = (eb > 0) ? (bg / eb) * (bg / eb) : 0.;

    char faixa[40];
    sprintf(faixa, "%.0f - %.0f",
            S->GetBinLowEdge(b), S->GetBinLowEdge(b) + S->GetBinWidth(b));

    printf("  %10s %11.4e %11.4e %7.0f", faixa, s, bg, nB);

    if (bg <= 0) printf(" %19s", "sem fundo MC");
    else         printf(" %12.1f%% +- %.1f%%",
                        100. * hist_pureza->GetBinContent(b),
                        100. * hist_pureza->GetBinError(b));

    if (hist_pureza_pos->GetBinContent(b) <= 0) printf(" %19s", "-");
    else printf(" %12.1f%% +- %.1f%%",
                100. * hist_pureza_pos->GetBinContent(b),
                100. * hist_pureza_pos->GetBinError(b));
    printf("\n");
  }

  // ---- numeros globais, integrando tudo de 5 a 40 GeV/c
  Int_t nb = hist_sinal->GetNbinsX();

  Double_t S_tot = hist_sinal->Integral(0, nb + 1);
  Double_t B_tot = hist_fundo->Integral(0, nb + 1);
  Double_t F_tot = hist_ruido[2]->Integral(0, nb + 1);

  Double_t S_pos = hist_sinal_pos->Integral(0, nb + 1);   // MEDIDO
  Double_t B_pos = hist_fundo_pos->Integral(0, nb + 1);

  printf("------------------------------------------\n");
  printf("  INTEGRADO de 5 a 40 GeV/c\n");
  printf("  %-22s %13s %13s\n", "", "so' isolamento", "+ massa");
  printf("  %-22s %13.4e %13.4e\n", "sinal (direto) [mb]", S_tot, S_pos);
  printf("  %-22s %13.4e %13.4e\n", "ruido total [mb]",    B_tot, B_pos);
  printf("  %-22s %13.4e %13.4e\n", "  pi0 + eta",
         hist_ruido[0]->Integral(0, nb + 1),
         hist_ruido_pos[0]->Integral(0, nb + 1));
  printf("  %-22s %13.4e %13.4e\n", "  outro hadron",
         hist_ruido[1]->Integral(0, nb + 1),
         hist_ruido_pos[1]->Integral(0, nb + 1));
  printf("  %-22s %13.4e %13.4e\n", "  fragmentacao",
         F_tot, hist_ruido_pos[2]->Integral(0, nb + 1));

  printf("------------------------------------------\n");
  if (S_tot + B_tot > 0)
  {
    printf("  PUREZA so' com isolamento   : %6.2f%%\n",
           100. * S_tot / (S_tot + B_tot));
    printf("  PUREZA com massa invariante : %6.2f%%\n",
           (S_pos + B_pos) > 0 ? 100. * S_pos / (S_pos + B_pos) : 0.);
    printf("  ganho                       : %6.2f pontos\n",
           ((S_pos + B_pos) > 0 ? 100. * S_pos / (S_pos + B_pos) : 0.)
           - 100. * S_tot / (S_tot + B_tot));
    printf("\n");
    printf("  O que o corte de massa fez:\n");
    printf("    ruido tirado  : %6.2f%%\n",
           B_tot ? 100. * (B_tot - B_pos) / B_tot : 0.);
    printf("    sinal perdido : %6.2f%%   (MEDIDO, nao suposto)\n",
           S_tot ? 100. * (S_tot - S_pos) / S_tot : 0.);
    printf("      A marcacao agora e' aplicada as DUAS amostras, com os\n");
    printf("      mesmos criterios. O sinal perdido e' o foton direto que\n");
    printf("      caiu numa das janelas de massa por acaso.\n");
    // O corte so' vale a pena se tirar MAIS ruido do que sinal. A razao
    // entre as duas fracoes diz quantas vezes mais: acima de 1 melhora a
    // pureza, abaixo de 1 piora.
    printf("    razao ruido tirado / sinal perdido : ");
    if (S_tot > S_pos && B_tot > 0)
      printf("%.1f\n", ((B_tot - B_pos) / B_tot) / ((S_tot - S_pos) / S_tot));
    else
      printf("sinal intacto, perda nula\n");
    printf("\n");
    printf("  PUREZA contando fragmentacao como sinal (so' isolamento):\n");
    printf("    %6.2f%%\n", 100. * (S_tot + F_tot) / (S_tot + B_tot));
  }

  // ---- quanta estatistica de Monte Carlo sustenta tudo isso
  {
    Double_t nB_cru = 0;
    for (Int_t o = 0; o < NORIGENS; o++) nB_cru += hist_ruido[o]->GetEntries();

    printf("------------------------------------------\n");
    printf("  ESTATISTICA DE MONTE CARLO\n");
    printf("    eventos no sinal : %.0f\n", hist_sinal->GetEntries());
    printf("    eventos no fundo : %.0f\n", nB_cru);

    if (nB_cru < 500)
    {
      printf("\n    AVISO: o fundo tem POUCOS eventos de Monte Carlo.\n");
      printf("    O valor integrado e' dominado pela regiao de pT baixo,\n");
      printf("    onde ha' estatistica, e la' ele e' confiavel. Mas a\n");
      printf("    CURVA bin a bin em pT alto nao e': um unico evento de\n");
      printf("    Monte Carlo carrega um peso grande, e os bins sem nenhum\n");
      printf("    evento de fundo ficam vazios no grafico em vez de\n");
      printf("    aparecer como pureza de 100%%, que seria falso.\n");
    }
  }
  printf("==========================================\n");

  // ================================================================
  //  6) SALVA
  // ================================================================
  TFile *arquivo = new TFile("pureza.root", "RECREATE");

  hist_sinal        ->SetDirectory(arquivo);
  hist_sinal_pos    ->SetDirectory(arquivo);
  hist_eiso_sinal   ->SetDirectory(arquivo);
  hist_fundo        ->SetDirectory(arquivo);
  hist_fundo_pos    ->SetDirectory(arquivo);
  hist_pureza       ->SetDirectory(arquivo);
  hist_pureza_pos   ->SetDirectory(arquivo);
  hist_pureza_prompt->SetDirectory(arquivo);

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    hist_ruido    [o]->SetDirectory(arquivo);
    hist_ruido_pos[o]->SetDirectory(arquivo);
    hist_eiso     [o]->SetDirectory(arquivo);
  }

  // Sem Close(): fechar o TFile destroi os histogramas que agora pertencem
  // a ele, e os canvases ficariam em branco. O Write() grava.
  arquivo->Write();
  printf("  Arquivo salvo: pureza.root\n");
  printf("  Figuras: pureza_espectros.png  pureza_espectros_norm.png\n");
  printf("           pureza_invariante.png  pureza_eiso.png\n");
  printf("           pureza_curva.png\n");
}
