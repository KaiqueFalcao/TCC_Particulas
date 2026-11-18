// Soma as 4 faixas EXCLUSIVAS de pThat, cada uma pesada pela SUA PROPRIA
// secao de choque, e redesenha todos os graficos da jatos_foton_v2.cc.
//
// RODAR
//   root -l -b -q soma_pthat.cc
//
// DE ONDE VEM O PESO: cada rodada da jatos_foton_v2.cc guarda, no bin 5 do
// histograma hist_sigma, o valor sigma/N -- a secao de choque daquela faixa
// de pThat dividida pelo numero de eventos aceitos. E' quantos milibarn vale
// um evento daquela rodada. Multiplicando o histograma por esse peso, o eixo
// vertical deixa de ser contagem e passa a ser secao de choque.
//
// POR QUE A SOMA ESTA' CERTA AQUI: as 4 faixas sao exclusivas (4-7, 7-10,
// 10-15 e acima de 15), ou seja, nenhum evento pertence a duas delas. Entao
//      h_total = soma_i ( sigma_i / N_i ) * h_i
// reconstroi o espectro de pThat > 4 inteiro, sem contar nada duas vezes, e
// com a estatistica reforcada na parte de pThat alto.
// NAO existe media de secoes de choque nesta conta: cada faixa entra com a
// sua propria sigma.
//
// Os histogramas da jatos_foton_v2.cc tem Sumw2 ligado, entao o Scale e o Add
// propagam as incertezas corretamente.
//
// UNIDADE DO EIXO VERTICAL: o padrao e' SECAO DE CHOQUE, em mb por bin. E' o
// que a soma pesada produz naturalmente, e e' a unidade em que o resultado se
// compara com dado experimental e com a foton_direto.cc. Quem quiser o eixo em
// numero de ocorrencias passa o numero de eventos equivalente:
//   root -l -b -q soma_pthat.cc                 // mb  (padrao)
//   root -l -b -q 'soma_pthat.cc(100000)'       // ocorrencias em 100 mil ev
//
// TRES MANEIRAS DE OLHAR O MESMO ESPECTRO, e todas aparecem nos PNGs:
//   1) secao de choque por bin, em mb            -- o padrao, ja' descrito
//   2) normalizado por foton, (1/N) dN/dX        -- area 1, so' a FORMA
//   3) secao de choque INVARIANTE,
//        (1/2pi) d2sigma/(dpT dy) em mb/(GeV/c)  -- a forma que se compara
//        direto com dado publicado, porque nao depende da largura do bin
//        nem do intervalo de rapidez usado
// A (2) e a (3) sao so' mudanca de escala da (1). NENHUMA DELAS entra na
// pureza: a pureza e' razao de secoes de choque, sempre a partir da (1).

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"

const Int_t NRODADAS = 4;
const Int_t NFAIXAS  = 6;
const Int_t NORIGENS = 3;   // 0 = pi0 + eta, 1 = outro hadron, 2 = fragmentacao

// Massas e janelas da marcacao, as mesmas da jatos_foton_v2.cc
const Double_t M_PI0      = 0.1349768;   // GeV
const Double_t M_ETA      = 0.547862;    // GeV
const Double_t JANELA_PI0 = 0.02;        // GeV
const Double_t JANELA_ETA = 0.02;        // GeV

// Corte de isolamento, so' para desenhar a linha vertical no grafico de E_iso
const Double_t PT_ISO_MAX = 2.0;         // GeV/c

// ============================================================================
//  PARAMETROS DA SECAO DE CHOQUE INVARIANTE  --  conferir com o orientador
// ============================================================================
// A forma invariante do espectro e'
//
//     (1/2pi) d2sigma/(dpT dy)  =  sigma_bin / ( 2pi * largura_do_bin * dy )
//
// em mb/(GeV/c). O dy e' o intervalo de rapidez em que os fotons foram
// selecionados: o corte e' |eta| < 0,5, e para foton (que nao tem massa)
// y = eta, entao dy = 0,5 - (-0,5) = 1,0.
//
// Ligando SP_USA_1_SOBRE_PT, entra tambem o fator 1/pT e a conta passa a ser
//
//     E d3sigma/dp3 = (1/(2pi pT)) d2sigma/(dpT dy)   em mb/(GeV/c)^2
//
// que e' a outra convencao usada na literatura. As duas descrevem a mesma
// medida; muda a potencia de pT no eixo. DEIXADO EM kFALSE ate' confirmar.
const Double_t SP_DY = 1.0;
const Bool_t   SP_USA_1_SOBRE_PT = kFALSE;

// As 4 rodadas, na ordem das faixas de pThat
const char *sp_sufixo[NRODADAS] = {"pthat04_07", "pthat07_10",
                                   "pthat10_15", "pthat15_inf"};

TFile    *sp_arquivo[NRODADAS];
Double_t  sp_peso   [NRODADAS];

// Quantos histogramas a soma ja' juntou. Contado em vez de escrito a mao:
// assim a mensagem do fim nunca fica desatualizada quando se acrescenta
// histograma novo.
Int_t     sp_nsomados = 0;

// Fator final aplicado a todos os somados, e o rotulo do eixo vertical.
// Em modo "ocorrencias" o fator converte mb de volta para contagem; em
// modo "mb" ele vale 1 e nada muda.
Double_t    sp_fator = 1.;
const char *sp_eixo  = "#sigma [mb] por bin";

// ---- histogramas SOMADOS. Mesmos nomes da jatos_foton_v2.cc, para a
// ---- correspondencia ficar obvia. Em escopo de arquivo, entao continuam
// ---- acessiveis no prompt do ROOT depois que a macro termina.
TH1F *hist_pt_jato    = 0;
TH1F *hist_n_part     = 0;
TH1F *hist_pt_foton   = 0;
TH2F *hist_eta_phi    = 0;
TH1F *hist_assimetria = 0;
TH1F *hist_z_gama     = 0;
TH1F *hist_soma_cone  = 0;

// Espectro de pT dos fotons ISOLADOS, um por origem (35 bins de 5 a 40
// GeV/c). Mesmo binning do hist_pt_direto_iso da foton_direto.cc.
TH1F *hist_pt_iso_pi0eta = 0;
TH1F *hist_pt_iso_frag   = 0;
TH1F *hist_pt_iso_outro  = 0;

// E_iso = soma de pT no cone, por origem, para todo foton energetico
TH1F *hist_eiso[NORIGENS];

// Massa invariante dos pares de fotons, nas duas janelas
TH1F *hist_massa_gg = 0;
TH1F *hist_massa_melhor    [NORIGENS];
TH1F *hist_massa_melhor_eta[NORIGENS];

// Espectro DEPOIS da marcacao: os marcados (pi0 ou eta) foram removidos
TH1F *hist_pt_pos_pi0eta = 0;
TH1F *hist_pt_pos_frag   = 0;
TH1F *hist_pt_pos_outro  = 0;

TH2F *hist_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F *hist_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F *hist_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};

TH1F *hist_z_iso   [NORIGENS][NFAIXAS];
TH1F *hist_z_naoiso[NORIGENS][NFAIXAS];

TCanvas *c1  = 0;
TCanvas *c2  = 0;
TCanvas *c3  = 0;
TCanvas *c4  = 0;
TCanvas *c5  = 0;
TCanvas *c6  = 0;
TCanvas *c7  = 0;
TCanvas *c8  = 0;
TCanvas *c9  = 0;
TCanvas *c10 = 0;
TCanvas *c11 = 0;
TCanvas *c12 = 0;
TCanvas *c13 = 0;   // espectro de pT dos isolados, por origem
TCanvas *c14 = 0;   // massa invariante
TCanvas *c15 = 0;   // gerais normalizados por foton
TCanvas *c16 = 0;   // E_iso por origem
TCanvas *c17 = 0;   // pT dos isolados, normalizado por foton
TCanvas *c18 = 0;   // pT dos isolados, secao de choque invariante

// Soma um histograma de UMA dimensao das 4 rodadas, cada uma multiplicada
// pelo peso da sua faixa. Devolve 0 se o histograma nao existir em algum
// dos arquivos.
//
// O clone recebe um nome temporario e so' depois e' rebatizado: se clonasse
// ja' com o nome final, ele nasceria dentro do TFile aberto, com nome
// repetido. O SetDirectory(0) solta o clone de qualquer arquivo, para que
// ele sobreviva ao Close la' no fim.
TH1F *soma_1d(const char *nome_hist)
{
  TH1F *total = 0;

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    TH1F *h = (TH1F*) sp_arquivo[i]->Get(nome_hist);

    if (!h)
    {
      printf("  ERRO: %s nao existe em %s\n", nome_hist, sp_sufixo[i]);
      return 0;
    }

    if (i == 0)
    {
      total = (TH1F*) h->Clone("sp_tmp_total");
      total->SetDirectory(0);
      total->SetName(nome_hist);
      total->Scale(sp_peso[0]);
    }
    else
    {
      TH1F *parcela = (TH1F*) h->Clone("sp_tmp_parcela");
      parcela->SetDirectory(0);
      parcela->Scale(sp_peso[i]);
      total->Add(parcela);
      delete parcela;
    }
  }

  // Depois do peso o conteudo e' secao de choque, nao contagem. O rotulo
  // do eixo tem que dizer isso, senao o grafico mente.
  total->Scale(sp_fator);
  total->SetYTitle(sp_eixo);
  sp_nsomados++;
  return total;
}

// Divide um espectro pela sua propria area, de modo que a area debaixo da
// curva passe a valer 1. O eixo deixa de ser secao de choque e passa a ser
// (1/N) dN/dX -- so' a FORMA da distribuicao.
//
// Scale(c, "width") divide cada bin pela largura dele e multiplica por c.
// Com c = 1/integral, sai exatamente a densidade de area 1.
//
// NAO serve para a pureza: razao de curvas normalizadas nao e' pureza.
TH1F *faz_normalizado(TH1F *h, const char *nome)
{
  TH1F *nm = (TH1F*) h->Clone(nome);
  nm->SetDirectory(0);                // copia de grafico: fora do .root

  // Integral(0, nbins+1) e' de proposito: inclui underflow e OVERFLOW.
  // A area tem que ser o numero TOTAL de fotons, nao so' os que couberam
  // no eixo. No E_iso, que vai de 0 a 20 GeV/c, 20% dos fotons de pi0+eta
  // e 16% dos de outro hadron tem E_iso maior que 20 e caem no overflow.
  // Com Integral() puro essas duas curvas subiam ~25% e ~18%, e o erro
  // NAO se cancela entre as curvas: o sinal e a fragmentacao quase nao tem
  // overflow, entao o grafico levantava so' o fundo.
  // O desenho continua de 0 a 20; o que muda e' que a area visivel agora
  // e' a FRACAO de fotons ali dentro, nao 1.
  Int_t nbx = nm->GetNbinsX();
  Double_t area = nm->Integral(0, nbx + 1);
  if (area > 0) nm->Scale(1. / area, "width");

  nm->SetStats(0);
  return nm;
}

// Converte um espectro de "sigma [mb] por bin" para SECAO DE CHOQUE
// INVARIANTE, dividindo cada bin por 2pi * largura * dy (e por pT, se a
// chave SP_USA_1_SOBRE_PT estiver ligada). O erro e' dividido pelo mesmo
// fator, que e' uma constante por bin.
TH1F *faz_invariante(TH1F *h, const char *nome)
{
  TH1F *inv = (TH1F*) h->Clone(nome);
  inv->SetDirectory(0);               // copia de grafico: fora do .root

  for (Int_t b = 1; b <= inv->GetNbinsX(); b++)
  {
    Double_t fator = 2. * TMath::Pi() * inv->GetBinWidth(b) * SP_DY;
    if (SP_USA_1_SOBRE_PT) fator *= inv->GetBinCenter(b);
    if (fator <= 0) continue;

    inv->SetBinContent(b, inv->GetBinContent(b) / fator);
    inv->SetBinError  (b, inv->GetBinError(b)   / fator);
  }

  inv->SetYTitle(SP_USA_1_SOBRE_PT
      ? "E d^{3}#sigma/dp^{3}  [mb/(GeV/c)^{2}]"
      : "(1/2#pi) d^{2}#sigma/(dp_{T}dy)  [mb/(GeV/c)]");
  inv->SetStats(0);
  return inv;
}

// Mesma coisa para os de DUAS dimensoes
TH2F *soma_2d(const char *nome_hist)
{
  TH2F *total = 0;

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    TH2F *h = (TH2F*) sp_arquivo[i]->Get(nome_hist);

    if (!h)
    {
      printf("  ERRO: %s nao existe em %s\n", nome_hist, sp_sufixo[i]);
      return 0;
    }

    if (i == 0)
    {
      total = (TH2F*) h->Clone("sp_tmp_total2");
      total->SetDirectory(0);
      total->SetName(nome_hist);
      total->Scale(sp_peso[0]);
    }
    else
    {
      TH2F *parcela = (TH2F*) h->Clone("sp_tmp_parcela2");
      parcela->SetDirectory(0);
      parcela->Scale(sp_peso[i]);
      total->Add(parcela);
      delete parcela;
    }
  }

  total->Scale(sp_fator);
  total->SetZTitle(sp_eixo);
  sp_nsomados++;
  return total;
}

// eventos_equiv = 0  : eixo em SECAO DE CHOQUE, mb por bin. E' o PADRAO.
// eventos_equiv > 0  : eixo em NUMERO DE OCORRENCIAS -- quantas haveria numa
//                      rodada inclusiva de pThat > 4 com esse tanto de
//                      eventos. Passe 100000 para o equivalente ao que cada
//                      faixa gerou.
void soma_pthat(Double_t eventos_equiv = 0.)
{
  const char *curto_origem[NORIGENS] = {"pi0eta", "outro", "frag"};
  const char *nome_origem [NORIGENS] = {"#pi^{0} + #eta", "outro hadron",
                                        "fragmentacao"};

  char nome[200], png[200];

  // ================================================================
  //  1) ABRE OS 4 ARQUIVOS E LE OS PESOS
  // ================================================================
  printf("\n============== RODADAS E PESOS ==============\n");
  printf("  %-13s %12s %12s %16s\n",
         "rodada", "sigma[mb]", "N", "peso[mb/ev]");

  Double_t sigma_total = 0;

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(nome, "jatos_foton_%s.root", sp_sufixo[i]);
    sp_arquivo[i] = TFile::Open(nome);

    if (!sp_arquivo[i] || sp_arquivo[i]->IsZombie())
    {
      printf("  ERRO: nao consegui abrir %s\n", nome);
      printf("  Rode as 4 producoes antes de somar:\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(4., 7.,43,100000,\"pthat04_07\")'\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(7.,10.,44,100000,\"pthat07_10\")'\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(10.,15.,45,100000,\"pthat10_15\")'\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(15.,-1.,46,100000,\"pthat15_inf\")'\n");
      return;
    }

    TH1D *info = (TH1D*) sp_arquivo[i]->Get("hist_sigma");

    if (!info)
    {
      printf("  ERRO: %s nao tem hist_sigma.\n", nome);
      printf("  Foi gerado por uma versao antiga da jatos_foton_v2.cc?\n");
      return;
    }

    sp_peso[i] = info->GetBinContent(5);   // bin 5 = sigma / n_aceitos
    sigma_total += info->GetBinContent(1);

    printf("  %-13s %12.4f %12.0f %16.6e\n", sp_sufixo[i],
           info->GetBinContent(1), info->GetBinContent(3), sp_peso[i]);
  }

  printf("---------------------------------------------\n");
  printf("  sigma total das 4 faixas : %.4f mb\n", sigma_total);
  printf("  Cada faixa entra com a SUA sigma. Nao ha' media de sigma.\n");

  // ---- escolhe a unidade do eixo vertical
  // A soma pesada nasce em mb. Para virar NUMERO DE OCORRENCIAS basta
  // dividir pela sigma total (dando ocorrencias por mb) e multiplicar pelo
  // tamanho da amostra equivalente. Com eventos_equiv = 100000, o resultado
  // e' "quantas ocorrencias haveria em 100 mil eventos de pThat > 4" -- o
  // mesmo eixo de uma rodada unica, mas com a estatistica das 4 faixas.
  if (eventos_equiv > 0 && sigma_total > 0)
  {
    sp_fator = eventos_equiv / sigma_total;
    sp_eixo  = "ocorrencias";
    printf("  EIXO EM OCORRENCIAS, equivalente a %.0f eventos de pThat > 4\n",
           eventos_equiv);
    printf("  (fator aplicado: %.4f ocorrencias por mb)\n", sp_fator);
  }
  else
  {
    sp_fator = 1.;
    sp_eixo  = "#sigma [mb] por bin";
    printf("  EIXO EM mb (milibarn), ou seja, secao de choque.\n");
  }

  // ================================================================
  //  2) SOMA TODOS OS HISTOGRAMAS
  //  O hist_sigma fica de fora: ele guarda informacao da rodada,
  //  nao e' um espectro para somar.
  // ================================================================
  hist_pt_jato    = soma_1d("hist_pt_jato");
  hist_n_part     = soma_1d("hist_n_part");
  hist_pt_foton   = soma_1d("hist_pt_foton");
  hist_assimetria = soma_1d("hist_assimetria");
  hist_z_gama     = soma_1d("hist_z_gama");
  hist_soma_cone  = soma_1d("hist_soma_cone");

  // espectro de pT dos isolados, por origem
  hist_pt_iso_pi0eta = soma_1d("hist_pt_iso_pi0eta");
  hist_pt_iso_frag   = soma_1d("hist_pt_iso_frag");
  hist_pt_iso_outro  = soma_1d("hist_pt_iso_outro");

  // E_iso por origem
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "hist_eiso_%s", curto_origem[o]);
    hist_eiso[o] = soma_1d(nome);
  }

  // massa invariante dos pares de fotons, nas duas janelas
  hist_massa_gg = soma_1d("hist_massa_gg");
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "hist_massa_melhor_%s", curto_origem[o]);
    hist_massa_melhor[o] = soma_1d(nome);

    sprintf(nome, "hist_massa_melhor_eta_%s", curto_origem[o]);
    hist_massa_melhor_eta[o] = soma_1d(nome);
  }

  // espectro depois da marcacao
  hist_pt_pos_pi0eta = soma_1d("hist_pt_pos_pi0eta");
  hist_pt_pos_frag   = soma_1d("hist_pt_pos_frag");
  hist_pt_pos_outro  = soma_1d("hist_pt_pos_outro");

  hist_eta_phi    = soma_2d("hist_eta_phi");

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    sprintf(nome, "hist_td_%d", f);  hist_td[f] = soma_2d(nome);
    sprintf(nome, "hist_Z_%d",  f);  hist_Z [f] = soma_2d(nome);
    sprintf(nome, "hist_AZ_%d", f);  hist_AZ[f] = soma_2d(nome);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      sprintf(nome, "hist_z_iso_%s_%d", curto_origem[o], f);
      hist_z_iso[o][f] = soma_1d(nome);

      sprintf(nome, "hist_z_naoiso_%s_%d", curto_origem[o], f);
      hist_z_naoiso[o][f] = soma_1d(nome);
    }
  }

  // Confere que nada falhou antes de desenhar
  if (!hist_pt_jato || !hist_eta_phi || !hist_z_iso[2][5])
  {
    printf("  Algum histograma nao foi somado. Parando aqui.\n");
    return;
  }

  // Os tres espectros por origem sao novos: se faltarem, o .root foi gerado
  // por uma versao da jatos_foton_v2.cc anterior a eles.
  if (!hist_pt_iso_pi0eta || !hist_pt_iso_frag || !hist_pt_iso_outro ||
      !hist_massa_gg       || !hist_massa_melhor[2] || !hist_pt_pos_pi0eta ||
      !hist_massa_melhor_eta[2] || !hist_eiso[2])
  {
    printf("  Faltam histogramas novos (hist_pt_iso_*, hist_eiso_* ou\n");
    printf("  hist_massa_*). O .root pode ser de uma versao anterior a\n");
    printf("  juncao pi0 + eta, que mudou o nome dos histogramas.\n");
    printf("  Rode a producao de novo:\n");
    printf("    ./producao.sh\n");
    return;
  }

  // ================================================================
  //  3) FECHA OS ARQUIVOS DE ENTRADA
  //  Pode fechar: todos os somados levaram SetDirectory(0) e estao
  //  na memoria, soltos de qualquer arquivo.
  // ================================================================
  for (Int_t i = 0; i < NRODADAS; i++) sp_arquivo[i]->Close();

  printf("  4 arquivos lidos e fechados. %d histogramas somados.\n",
         sp_nsomados);

  // ================================================================
  //  4) ARQUIVO DE SAIDA
  //  Criado agora, DEPOIS da soma. Os somados sao entregues a ele
  //  com SetDirectory, para que o Write() os grave.
  // ================================================================
  TFile *arquivo = new TFile("soma_pthat.root", "RECREATE");

  hist_pt_jato   ->SetDirectory(arquivo);
  hist_n_part    ->SetDirectory(arquivo);
  hist_pt_foton  ->SetDirectory(arquivo);
  hist_assimetria->SetDirectory(arquivo);
  hist_z_gama    ->SetDirectory(arquivo);
  hist_soma_cone ->SetDirectory(arquivo);
  hist_pt_iso_pi0eta->SetDirectory(arquivo);
  hist_pt_iso_frag  ->SetDirectory(arquivo);
  hist_pt_iso_outro ->SetDirectory(arquivo);
  hist_massa_gg     ->SetDirectory(arquivo);
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    hist_eiso[o]                ->SetDirectory(arquivo);
    hist_massa_melhor[o]        ->SetDirectory(arquivo);
    hist_massa_melhor_eta[o]    ->SetDirectory(arquivo);
  }
  hist_pt_pos_pi0eta->SetDirectory(arquivo);
  hist_pt_pos_frag  ->SetDirectory(arquivo);
  hist_pt_pos_outro ->SetDirectory(arquivo);
  hist_eta_phi   ->SetDirectory(arquivo);

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    hist_td[f]->SetDirectory(arquivo);
    hist_Z [f]->SetDirectory(arquivo);
    hist_AZ[f]->SetDirectory(arquivo);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      hist_z_iso   [o][f]->SetDirectory(arquivo);
      hist_z_naoiso[o][f]->SetDirectory(arquivo);
    }
  }

  // ================================================================
  //  5) CANVASES -- os mesmos da jatos_foton_v2.cc, com prefixo soma_
  // ================================================================

  // ---- c1: histogramas gerais
  c1 = new TCanvas("c1", "SOMA -- jatos e fotons, gerais", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato ->Draw("HIST E");
  c1->cd(2); hist_n_part  ->Draw("HIST E");
  c1->cd(3); hist_pt_foton->Draw("HIST E");
  c1->cd(4); gPad->SetLogy(); hist_soma_cone->Draw("HIST E");
  c1->Update();
  c1->SaveAs("soma_gerais.png");

  // ---- c2: mapa eta-phi
  c2 = new TCanvas("c2", "SOMA -- mapa eta-phi dos fotons", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();
  c2->SaveAs("soma_mapa_eta_phi.png");

  // ---- c3: assimetria e foton lider
  c3 = new TCanvas("c3", "SOMA -- fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw("HIST E");
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw("HIST E");
  c3->Update();
  c3->SaveAs("soma_assimetria.png");

  // ---- c4, c5, c6: lider x nao-lider, os 3 niveis de corte
  c4 = new TCanvas("c4", "SOMA -- lider x nao-lider (todos)", 1500, 1000);
  c4->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c4->cd(f + 1); hist_td[f]->Draw("COLZ"); }
  c4->Update();
  c4->SaveAs("soma_lider_todos.png");

  c5 = new TCanvas("c5", "SOMA -- lider x nao-lider (z alto)", 1500, 1000);
  c5->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c5->cd(f + 1); hist_Z[f]->Draw("COLZ"); }
  c5->Update();
  c5->SaveAs("soma_lider_z.png");

  c6 = new TCanvas("c6", "SOMA -- lider x nao-lider (z e A altos)", 1500, 1000);
  c6->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c6->cd(f + 1); hist_AZ[f]->Draw("COLZ"); }
  c6->Update();
  c6->SaveAs("soma_lider_z_e_assimetria.png");

  // ---- c7, c8, c9: z dos ISOLADOS, uma origem por canvas
  TCanvas *canvas_iso[NORIGENS];
  const char *arq_iso[NORIGENS] = {"soma_z_iso_pi0eta.png",
                                   "soma_z_iso_outro.png",
                                   "soma_z_iso_frag.png"};
  const char *tit_iso[NORIGENS] = {"SOMA -- z: isolado, de pi0 + eta",
                                   "SOMA -- z: isolado, de outro hadron",
                                   "SOMA -- z: isolado, de fragmentacao"};

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "c%d", 7 + o);
    canvas_iso[o] = new TCanvas(nome, tit_iso[o], 1500, 1000);
    canvas_iso[o]->Divide(3, 2);
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      canvas_iso[o]->cd(f + 1);
      gPad->SetLeftMargin(0.16);          // senao os numeros saem cortados
      hist_z_iso[o][f]->Draw("HIST E");
    }
    canvas_iso[o]->Update();
    canvas_iso[o]->SaveAs(arq_iso[o]);
  }
  c7 = canvas_iso[0];  c8 = canvas_iso[1];  c9 = canvas_iso[2];

  // ---- c10, c11, c12: z dos NAO ISOLADOS, para comparar
  TCanvas *canvas_nao[NORIGENS];
  const char *arq_nao[NORIGENS] = {"soma_z_naoiso_pi0eta.png",
                                   "soma_z_naoiso_outro.png",
                                   "soma_z_naoiso_frag.png"};
  const char *tit_nao[NORIGENS] = {"SOMA -- z: nao isolado, de pi0 + eta",
                                   "SOMA -- z: nao isolado, de outro hadron",
                                   "SOMA -- z: nao isolado, de fragmentacao"};

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "c%d", 10 + o);
    canvas_nao[o] = new TCanvas(nome, tit_nao[o], 1500, 1000);
    canvas_nao[o]->Divide(3, 2);
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      canvas_nao[o]->cd(f + 1);
      gPad->SetLeftMargin(0.16);
      hist_z_naoiso[o][f]->Draw("HIST E");
    }
    canvas_nao[o]->Update();
    canvas_nao[o]->SaveAs(arq_nao[o]);
  }
  c10 = canvas_nao[0];  c11 = canvas_nao[1];  c12 = canvas_nao[2];

  // ---- c13: espectro de pT dos fotons ISOLADOS, um painel por origem.
  // Este e' o RUIDO separado por origem, no mesmo binning do sinal da
  // foton_direto.cc -- da' para por os quatro no mesmo grafico depois.
  TH1F *espectro[NORIGENS] = {hist_pt_iso_pi0eta, hist_pt_iso_frag,
                              hist_pt_iso_outro};
  const char *tit_esp[NORIGENS] = {"#gamma isolado de #pi^{0} + #eta",
                                   "#gamma isolado de fragmentacao",
                                   "#gamma isolado de outro hadron"};
  Int_t cor_esp[NORIGENS] = {kBlue, kRed, kGreen + 2};

  c13 = new TCanvas("c13", "SOMA -- pT dos isolados por origem", 1500, 500);
  c13->Divide(3, 1);

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    c13->cd(o + 1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.16);        // senao os expoentes saem cortados

    espectro[o]->SetTitle(tit_esp[o]);
    espectro[o]->SetXTitle("p_{T}^{#gamma} [GeV/c]");
    espectro[o]->SetYTitle(sp_eixo);
    espectro[o]->SetLineColor(cor_esp[o]);
    espectro[o]->SetLineWidth(2);
    espectro[o]->SetStats(0);
    espectro[o]->Draw("HIST E");
  }
  c13->Update();
  c13->SaveAs("soma_pt_iso_por_origem.png");

  // ---- c14: MASSA INVARIANTE, ja' somada e pesada por secao de choque.
  // Painel 1: todos os pares, alcance inteiro -- o pi0 e o eta aparecem.
  // Painel 2: zoom na regiao do pi0.
  // Painel 3: melhor par de cada foton isolado NA JANELA DO PI0, por origem.
  // Painel 4: o mesmo NA JANELA DO ETA.
  // Os paineis 3 e 4 sao o retrato da marcacao: mostram o que cada uma das
  // duas janelas pega em cada origem -- e, na curva de fragmentacao, o que
  // cada uma pega por engano.
  c14 = new TCanvas("c14", "SOMA -- massa invariante", 1400, 1000);
  c14->Divide(2, 2);

  c14->cd(1);
  gPad->SetLogy();  gPad->SetLeftMargin(0.16);
  hist_massa_gg->SetStats(0);
  hist_massa_gg->SetLineColor(kBlack);
  hist_massa_gg->SetLineWidth(2);
  hist_massa_gg->SetYTitle(sp_eixo);
  hist_massa_gg->Draw("HIST E");

  // zoom numa COPIA, para nao mexer no eixo do que vai para o .root
  c14->cd(2);
  gPad->SetLeftMargin(0.16);
  TH1F *zoom = (TH1F*) hist_massa_gg->Clone("sp_massa_zoom");
  zoom->SetDirectory(0);
  zoom->SetTitle("zoom no #pi^{0}");
  zoom->GetXaxis()->SetRangeUser(0., 0.3);
  zoom->SetLineColor(kBlue);
  zoom->Draw("HIST E");

  // Os dois paineis de baixo tem a mesma estrutura, so' muda o conjunto de
  // histogramas e a janela mostrada. Por isso vao num laco de dois passos.
  {
    Int_t cor_m[NORIGENS] = {kBlue, kGreen + 2, kRed};

    for (Int_t jan = 0; jan < 2; jan++)      // 0 = janela do pi0, 1 = do eta
    {
      TH1F **conj = (jan == 0) ? hist_massa_melhor : hist_massa_melhor_eta;
      Double_t x_lo = (jan == 0) ? 0.00 : 0.35;
      Double_t x_hi = (jan == 0) ? 0.30 : 0.75;

      c14->cd(3 + jan);
      gPad->SetLogy();  gPad->SetLeftMargin(0.16);

      Double_t topo_m = 0;
      for (Int_t o = 0; o < NORIGENS; o++)
        if (conj[o]->GetMaximum() > topo_m) topo_m = conj[o]->GetMaximum();

      TLegend *leg_m = new TLegend(0.45, 0.70, 0.89, 0.89);
      leg_m->SetTextSize(0.035);

      for (Int_t o = 0; o < NORIGENS; o++)
      {
        conj[o]->SetTitle(jan == 0 ? "melhor par na janela do #pi^{0}"
                                   : "melhor par na janela do #eta");
        conj[o]->SetStats(0);
        conj[o]->SetLineColor(cor_m[o]);
        conj[o]->SetLineWidth(2);
        conj[o]->SetMaximum(topo_m * 2.);
        conj[o]->SetYTitle(sp_eixo);
        conj[o]->GetXaxis()->SetRangeUser(x_lo, x_hi);
        conj[o]->Draw(o == 0 ? "HIST E" : "HIST E SAME");
        leg_m->AddEntry(conj[o], nome_origem[o], "l");
      }
      leg_m->Draw();
    }
  }
  c14->Update();
  c14->SaveAs("soma_massa_invariante.png");

  // ---- c15: os MESMOS graficos do soma_gerais.png, mas normalizados por
  // foton (ou por jato): cada curva dividida pela sua area, que passa a
  // valer 1. O eixo deixa de ser secao de choque e passa a ser (1/N) dN/dX.
  // Serve para comparar FORMAS; os tamanhos relativos ficam no soma_gerais.
  c15 = new TCanvas("c15", "SOMA -- gerais, normalizados", 1200, 1000);
  c15->Divide(2, 2);
  {
    TH1F *orig[4] = {hist_pt_jato, hist_n_part, hist_pt_foton, hist_soma_cone};
    const char *eixo_y[4] = {"(1/N) dN/dp_{T}^{jato}  [c/GeV]",
                             "(1/N) dN/dN_{const}",
                             "(1/N) dN/dp_{T}^{#gamma}  [c/GeV]",
                             "(1/N) dN/dE_{iso}  [c/GeV]"};

    for (Int_t i = 0; i < 4; i++)
    {
      c15->cd(i + 1);
      gPad->SetLeftMargin(0.16);
      if (i == 3) gPad->SetLogy();

      sprintf(nome, "sp_norm_%d", i);
      TH1F *nm = faz_normalizado(orig[i], nome);
      nm->SetYTitle(eixo_y[i]);
      nm->SetLineWidth(2);
      nm->Draw("HIST E");
    }
  }
  c15->Update();
  c15->SaveAs("soma_gerais_norm.png");

  // ---- c16: E_iso por origem, normalizado por foton (area 1 em cada).
  // A linha vertical tracejada e' o corte de isolamento: o que esta' a
  // esquerda e' aceito. A distancia entre as curvas nessa regiao e' o
  // motivo de o isolamento funcionar.
  c16 = new TCanvas("c16", "SOMA -- E_iso por origem", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  {
    Int_t cor_e[NORIGENS] = {kBlue, kGreen + 2, kRed};
    TH1F *norm_e[NORIGENS];
    Double_t topo_e = 0, piso_e = 1e30;

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      sprintf(nome, "sp_eiso_norm_%s", curto_origem[o]);
      norm_e[o] = faz_normalizado(hist_eiso[o], nome);

      if (norm_e[o]->GetMaximum() > topo_e) topo_e = norm_e[o]->GetMaximum();
      for (Int_t b = 1; b <= norm_e[o]->GetNbinsX(); b++)
      {
        Double_t v = norm_e[o]->GetBinContent(b);
        if (v > 0 && v < piso_e) piso_e = v;
      }
    }
    if (piso_e > 1e29) piso_e = 1e-6;

    TLegend *leg_e = new TLegend(0.48, 0.68, 0.89, 0.89);
    leg_e->SetTextSize(0.033);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      norm_e[o]->SetTitle("isolamento por origem, normalizado por #gamma");
      norm_e[o]->SetXTitle("E_{iso} = #Sigma pT^{cone} [GeV/c]");
      norm_e[o]->SetYTitle("(1/N_{#gamma}) dN/dE_{iso}  [c/GeV]");
      norm_e[o]->SetLineColor(cor_e[o]);
      norm_e[o]->SetLineWidth(2);
      norm_e[o]->SetMinimum(piso_e * 0.5);
      norm_e[o]->SetMaximum(topo_e * 3.);
      norm_e[o]->Draw(o == 0 ? "HIST E" : "HIST E SAME");
      leg_e->AddEntry(norm_e[o], nome_origem[o], "l");
    }

    TLine *lin_e = new TLine(PT_ISO_MAX, piso_e * 0.5,
                             PT_ISO_MAX, topo_e * 3.);
    lin_e->SetLineColor(kBlack);
    lin_e->SetLineStyle(2);
    lin_e->SetLineWidth(2);
    lin_e->Draw();
    leg_e->AddEntry(lin_e, "corte de isolamento", "l");
    leg_e->Draw();
  }
  c16->Update();
  c16->SaveAs("soma_eiso.png");

  // ---- c17: espectro de pT dos isolados, normalizado por foton.
  // Mesma informacao do soma_pt_iso_por_origem.png, outra escala: aqui as
  // tres curvas tem area 1 e da' para ver que o espectro do fundo de
  // decaimento cai mais rapido em pT do que o de fragmentacao.
  c17 = new TCanvas("c17", "SOMA -- pT dos isolados, normalizado", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  {
    TH1F *bruto[NORIGENS] = {hist_pt_iso_pi0eta, hist_pt_iso_outro,
                             hist_pt_iso_frag};
    Int_t cor_n[NORIGENS] = {kBlue, kGreen + 2, kRed};

    TLegend *leg_n = new TLegend(0.48, 0.72, 0.89, 0.89);
    leg_n->SetTextSize(0.033);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      sprintf(nome, "sp_pt_iso_norm_%s", curto_origem[o]);
      TH1F *nm = faz_normalizado(bruto[o], nome);

      nm->SetTitle("#gamma isolado por origem, normalizado por #gamma");
      nm->SetXTitle("p_{T}^{#gamma} [GeV/c]");
      nm->SetYTitle("(1/N_{#gamma}) dN/dp_{T}  [c/GeV]");
      nm->SetLineColor(cor_n[o]);
      nm->SetLineWidth(2);
      nm->Draw(o == 0 ? "HIST E" : "HIST E SAME");
      leg_n->AddEntry(nm, nome_origem[o], "l");
    }
    leg_n->Draw();
  }
  c17->Update();
  c17->SaveAs("soma_pt_iso_norm.png");

  // ---- c18: espectro de pT dos isolados em SECAO DE CHOQUE INVARIANTE,
  // antes e depois da marcacao por massa. Preto: a soma das tres origens.
  // Tracejado: a mesma soma depois de jogar fora os marcados. A diferenca
  // entre as duas e' o que a marcacao conseguiu remover.
  c18 = new TCanvas("c18", "SOMA -- pT invariante", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.17);
  {
    // soma das tres origens, antes e depois
    TH1F *antes = (TH1F*) hist_pt_iso_pi0eta->Clone("sp_fundo_antes");
    antes->SetDirectory(0);
    antes->Add(hist_pt_iso_outro);
    antes->Add(hist_pt_iso_frag);

    TH1F *depois = (TH1F*) hist_pt_pos_pi0eta->Clone("sp_fundo_depois");
    depois->SetDirectory(0);
    depois->Add(hist_pt_pos_outro);
    depois->Add(hist_pt_pos_frag);

    TH1F *inv_a = faz_invariante(antes,  "sp_inv_antes");
    TH1F *inv_d = faz_invariante(depois, "sp_inv_depois");

    inv_a->SetTitle("ruido isolado -- secao de choque invariante");
    inv_a->SetXTitle("p_{T}^{#gamma} [GeV/c]");
    inv_a->SetLineColor(kBlack);
    inv_a->SetLineWidth(3);
    inv_a->Draw("HIST E");

    inv_d->SetLineColor(kAzure + 2);
    inv_d->SetLineWidth(3);
    inv_d->SetLineStyle(2);
    inv_d->Draw("HIST E SAME");

    TLegend *leg_i = new TLegend(0.42, 0.72, 0.89, 0.89);
    leg_i->SetTextSize(0.032);
    leg_i->AddEntry(inv_a, "so' isolamento", "l");
    leg_i->AddEntry(inv_d, "isolamento + massa m_{#gamma#gamma}", "l");
    sprintf(nome, "dy = %.1f%s", SP_DY,
            SP_USA_1_SOBRE_PT ? ", com 1/p_{T}" : ", sem 1/p_{T}");
    leg_i->SetHeader(nome);
    leg_i->Draw();

    delete antes;  delete depois;
  }
  c18->Update();
  c18->SaveAs("soma_pt_iso_invariante.png");

  // ================================================================
  //  RESUMO
  // ================================================================
  // Os somados ja' foram multiplicados por sp_fator. Para voltar a mb basta
  // desfazer essa multiplicacao -- assim a tabela sai nas DUAS unidades e
  // nao fica dependendo de em que modo a macro foi chamada.
  Bool_t outra_unidade = (sp_fator != 1.);

  printf("=============================================\n");
  printf("  ESPECTRO SOMADO (integral)\n");
  printf("    pT dos jatos    : %.6e mb", hist_pt_jato->Integral() / sp_fator);
  if (outra_unidade) printf("   = %.6g %s", hist_pt_jato->Integral(), sp_eixo);
  printf("\n");
  printf("    pT dos fotons   : %.6e mb", hist_pt_foton->Integral() / sp_fator);
  if (outra_unidade) printf("   = %.6g %s", hist_pt_foton->Integral(), sp_eixo);
  printf("\n");

  printf("---------------------------------------------\n");
  printf("  FOTONS ENERGETICOS SOMADOS, por origem\n");

  if (outra_unidade)
    printf("  %-8s %14s %14s %13s %13s\n", "origem",
           "iso [mb]", "nao-iso [mb]", "iso [eixo]", "nao-iso[eixo]");
  else
    printf("  %-8s %14s %14s\n", "origem", "iso [mb]", "nao-iso [mb]");

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    Double_t s_iso = 0, s_nao = 0;
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      // Integral(0, nbins+1) e' de proposito: inclui o underflow e o
      // OVERFLOW. O histograma de z vai de 0 a 1, e um foton de
      // fragmentacao isolado e' praticamente o jato inteiro, ou seja
      // z ~ 1, que cai no overflow. Com Integral() puro esses fotons
      // ficavam de fora e a fracao isolada da fragmentacao saia baixa
      // (67% em vez de 73%). Nao e' corte nenhum: e' so' nao jogar
      // fora a ultima caixa.
      Int_t nb = hist_z_iso[o][f]->GetNbinsX();
      s_iso += hist_z_iso   [o][f]->Integral(0, nb + 1);
      s_nao += hist_z_naoiso[o][f]->Integral(0, nb + 1);
    }

    if (outra_unidade)
      printf("  %-8s %14.6e %14.6e %13.4f %13.4f\n", curto_origem[o],
             s_iso / sp_fator, s_nao / sp_fator, s_iso, s_nao);
    else
      printf("  %-8s %14.6e %14.6e\n", curto_origem[o],
             s_iso / sp_fator, s_nao / sp_fator);
  }

  // ---- MARCACAO POR MASSA INVARIANTE, pesada por secao de choque.
  // Esta e' a versao que vale: os contadores crus que cada rodada
  // imprime sao por faixa de pThat e nao podem ser somados direto.
  printf("---------------------------------------------\n");
  printf("  MARCACAO POR MASSA INVARIANTE (duas janelas)\n");
  printf("  Marcado = existe parceiro visivel com\n");
  printf("      |m - %.4f| < %.0f MeV   (pi0)   OU\n", M_PI0, JANELA_PI0 * 1000.);
  printf("      |m - %.4f| < %.0f MeV   (eta)\n",      M_ETA, JANELA_ETA * 1000.);
  // ATENCAO ao denominador. O hist_massa_melhor so' e' preenchido quando
  // existe ALGUM parceiro visivel. O foton cujo irmao ficou abaixo do
  // limiar nao entra nele -- e e' justamente esse o fundo IRREDUTIVEL.
  // Usar o hist_massa_melhor como denominador daria a eficiencia "entre
  // os que tem parceiro", que e' otimista e nao responde nada. O
  // denominador certo e' o total de isolados, do hist_pt_iso_*.
  printf("  %-8s %12s %12s %11s %11s %9s\n", "origem",
         "isolado[mb]", "c/ parc[mb]", "pi0 [mb]", "eta [mb]", "s/ parc%");

  {
    // na ordem de curto_origem: pi0 + eta, outro hadron, fragmentacao
    TH1F *iso_org[NORIGENS] = {hist_pt_iso_pi0eta, hist_pt_iso_outro,
                               hist_pt_iso_frag};

    Int_t p_lo = hist_massa_melhor[0]->FindBin(M_PI0 - JANELA_PI0);
    Int_t p_hi = hist_massa_melhor[0]->FindBin(M_PI0 + JANELA_PI0);
    Int_t e_lo = hist_massa_melhor_eta[0]->FindBin(M_ETA - JANELA_ETA);
    Int_t e_hi = hist_massa_melhor_eta[0]->FindBin(M_ETA + JANELA_ETA);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      // Integral(0, n+1) inclui o overflow: o foton de pT acima de 40
      // GeV/c existe e tem que entrar na conta.
      Int_t nb = iso_org[o]->GetNbinsX();
      Double_t tot  = iso_org[o]->Integral(0, nb + 1) / sp_fator;

      Double_t comp  = hist_massa_melhor[o]->Integral() / sp_fator;
      Double_t m_pi0 = hist_massa_melhor    [o]->Integral(p_lo, p_hi) / sp_fator;
      Double_t m_eta = hist_massa_melhor_eta[o]->Integral(e_lo, e_hi) / sp_fator;

      printf("  %-8s %12.4e %12.4e %11.4e %11.4e %8.2f%%\n",
             curto_origem[o], tot, comp, m_pi0, m_eta,
             tot ? 100. * (tot - comp) / tot : 0.);
    }

    printf("  (pi0 e eta sao as duas janelas separadas. A soma das duas e'\n");
    printf("   um LIMITE DE CIMA do total removido, porque um mesmo foton\n");
    printf("   pode ter um parceiro em cada janela e seria contado duas\n");
    printf("   vezes. O total que vale e' o da tabela seguinte, que vem do\n");
    printf("   que foi de fato removido.)\n");
    printf("  (s/ parc%% = sem parceiro visivel nenhum: o fundo IRREDUTIVEL,\n");
    printf("   o foton cujo irmao ficou abaixo do limiar.)\n");
  }

  // ---- confere a remocao contra os espectros de fato removidos.
  // A conta acima e' feita a partir do histograma de massa; esta e' feita
  // a partir dos hist_pt_pos_*, que o jatos_foton_v2.cc preencheu so' com
  // os fotons que SOBREVIVERAM. As duas tem que bater -- se nao baterem,
  // ha' desacordo entre o que foi contado e o que foi removido.
  printf("  ------------------------------------------\n");
  printf("  ESPECTRO REMOVIDO DE FATO (hist_pt_pos_*)\n");
  printf("  Esta e' a conta que vale: o jatos_foton_v2.cc so' preencheu o\n");
  printf("  hist_pt_pos_* com os fotons que SOBREVIVERAM as duas janelas.\n");
  printf("  %-8s %13s %13s %10s\n",
         "origem", "antes[mb]", "depois[mb]", "removido%");
  {
    TH1F *ant[NORIGENS] = {hist_pt_iso_pi0eta, hist_pt_iso_outro,
                           hist_pt_iso_frag};
    TH1F *dep[NORIGENS] = {hist_pt_pos_pi0eta, hist_pt_pos_outro,
                           hist_pt_pos_frag};

    Double_t antes_tot = 0, depois_tot = 0, frag_a = 0, frag_d = 0;

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      Int_t nb = ant[o]->GetNbinsX();
      Double_t a = ant[o]->Integral(0, nb + 1) / sp_fator;
      Double_t d = dep[o]->Integral(0, nb + 1) / sp_fator;

      printf("  %-8s %13.4e %13.4e %9.2f%%\n", curto_origem[o], a, d,
             a ? 100. * (a - d) / a : 0.);

      antes_tot += a;  depois_tot += d;
      if (o == 2) { frag_a = a; frag_d = d; }
    }

    printf("  ------------------------------------------\n");
    printf("  amostra isolada ANTES  : %.6e mb\n", antes_tot);
    printf("  amostra isolada DEPOIS : %.6e mb\n", depois_tot);
    printf("  rejeitado pela marcacao: %6.2f%%\n",
           antes_tot ? 100. * (antes_tot - depois_tot) / antes_tot : 0.);
    printf("  fragmentacao ANTES     : %6.2f%%\n",
           antes_tot  ? 100. * frag_a / antes_tot  : 0.);
    printf("  fragmentacao DEPOIS    : %6.2f%%\n",
           depois_tot ? 100. * frag_d / depois_tot : 0.);
  }

  // ---- E_ISO por origem: media e eficiencia do isolamento.
  // A fracao abaixo do corte e' a eficiencia do isolamento para aquela
  // origem. O corte de 2,0 GeV/c cai numa borda de bin, entao a integral
  // ate' o bin anterior e' exatamente o que e' aceito.
  printf("---------------------------------------------\n");
  printf("  E_iso POR ORIGEM (soma de pT no cone)\n");
  printf("  %-8s %13s %10s %13s %11s\n", "origem",
         "energet.[mb]", "<E_iso>", "aceito[mb]", "eficiencia");
  {
    Int_t b_corte = hist_eiso[0]->FindBin(PT_ISO_MAX) - 1;   // ultimo aceito

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      Int_t nbe = hist_eiso[o]->GetNbinsX();
      Double_t tudo   = hist_eiso[o]->Integral(0, nbe + 1) / sp_fator;
      Double_t aceito = hist_eiso[o]->Integral(1, b_corte)  / sp_fator;

      printf("  %-8s %13.4e %10.3f %13.4e %10.2f%%\n", curto_origem[o],
             tudo, hist_eiso[o]->GetMean(), aceito,
             tudo ? 100. * aceito / tudo : 0.);
    }
  }

  // ---- espectro de pT dos isolados, por origem (os histogramas novos).
  // Esta integral conta TODO foton isolado, inclusive os sem jato, entao e'
  // naturalmente maior que a coluna "iso" da tabela de cima, que vem dos
  // histogramas de z e por isso exige um jato associado.
  printf("---------------------------------------------\n");
  printf("  ESPECTRO DE pT DOS ISOLADOS, 5 a 40 GeV/c\n");
  printf("  (conta tambem os que nao cairam em jato)\n");

  TH1F *esp_res[NORIGENS] = {hist_pt_iso_pi0eta, hist_pt_iso_frag,
                             hist_pt_iso_outro};
  const char *nome_res[NORIGENS] = {"pi0eta", "frag", "outro"};

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    printf("  %-8s %14.6e mb", nome_res[o],
           esp_res[o]->Integral() / sp_fator);
    if (outra_unidade) printf("   = %.4f %s", esp_res[o]->Integral(), sp_eixo);
    printf("   (%.0f entradas brutas)\n", esp_res[o]->GetEntries());
  }

  printf("=============================================\n");

  // Sem Close(): fechar o TFile de saida destroi os histogramas que agora
  // pertencem a ele, e os canvases ficariam em branco. O Write() grava.
  arquivo->Write();
  printf("  Arquivo salvo: soma_pthat.root\n");
  printf("  18 PNGs salvos com prefixo soma_\n");
  printf("    em secao de choque : soma_gerais, soma_pt_iso_por_origem, ...\n");
  printf("    normalizados       : soma_gerais_norm, soma_pt_iso_norm,\n");
  printf("                         soma_eiso\n");
  printf("    invariante         : soma_pt_iso_invariante\n");
}
