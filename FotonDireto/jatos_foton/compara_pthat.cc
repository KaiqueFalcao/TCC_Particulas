// Combina as rodadas fatiadas em pThat.
//
// Le os 4 arquivos jatos_foton_pthat*.root, pega a secao de choque de cada um
// dentro do hist_sigma, e produz:
//   cs1) as 4 rodadas normalizadas e SOBREPOSTAS   (teste de consistencia)
//   cs2) a SOMA ponderada das 4
//   cs3) as 4 secoes de choque guardadas num histograma
//
// O PESO de cada rodada e'  w_i = sigma_i / N_i , em mb por evento. Ao
// multiplicar o histograma por w_i, o eixo vertical deixa de ser contagem e
// passa a ser secao de choque: quantos mb por bin de pT.
//
// NORMALIZAR PELA MAIOR SECAO DE CHOQUE: dividir todos os pesos por w da
// rodada de maior sigma (a de pThat > 4) deixa essa rodada com peso 1 e as
// outras em escala relativa. E' so' um fator global: nao muda a forma de
// nenhum espectro nem a razao entre eles. Ligue com normalizar_pela_maior.
//
// POR QUE A SOBREPOSICAO (cs1) IMPORTA: as 4 amostras sao ENCAIXADAS, porque
// pTHatMin = 4 gera tambem os eventos de pThat > 15. Normalizada, cada rodada
// da' o espectro CORRETO na regiao que ela cobre. Entao, onde todas valem (pT
// alto), as 4 curvas tem que cair UMA SOBRE A OUTRA. Se elas coincidem, somar
// as 4 conta aquela regiao 4 vezes -- e o cs2 mostra o tamanho do estrago.

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

const Int_t NRODADAS = 4;

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina.
TH1F   *esc_pt_jato [NRODADAS];   // pT do jato de cada rodada, ja' pesado
TH1F   *esc_pt_foton[NRODADAS];   // pT do foton de cada rodada, ja' pesado
TH1F   *soma_pt_jato  = 0;        // a soma ponderada das 4
TH1F   *soma_pt_foton = 0;
TH1D   *hist_sigmas   = 0;        // as 4 secoes de choque num histograma

TCanvas *cs1 = 0;
TCanvas *cs2 = 0;
TCanvas *cs3 = 0;

void compara_pthat(Bool_t normalizar_pela_maior = kTRUE)
{
  const char *sufixo[NRODADAS] = {"pthat04", "pthat07", "pthat10", "pthat15"};
  Int_t cor[NRODADAS] = {kBlack, kBlue, kGreen + 2, kRed};

  char nome[200], titulo[300], rotulo[60];

  // ---- arquivo de saida, criado ANTES dos clones para que eles pertencam
  // ---- a ele e entrem no Write()
  TFile *arquivo = new TFile("jatos_foton_combinado.root", "RECREATE");

  // ---- histograma com as 4 secoes de choque, um bin por rodada
  hist_sigmas = new TH1D("hist_sigmas",
      "secao de choque de cada rodada;;#sigma [mb]", NRODADAS, 0, NRODADAS);

  // ---- le as 4 rodadas
  TFile    *entrada[NRODADAS];
  Double_t  sigma[NRODADAS], sigma_err[NRODADAS], peso[NRODADAS];
  Double_t  n_aceitos[NRODADAS];

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(nome, "jatos_foton_%s.root", sufixo[i]);
    entrada[i] = TFile::Open(nome);

    if (!entrada[i] || entrada[i]->IsZombie())
    {
      printf("  ERRO: nao consegui abrir %s\n", nome);
      printf("  Rode as 4 producoes antes de combinar.\n");
      return;
    }

    TH1D *info = (TH1D*) entrada[i]->Get("hist_sigma");
    if (!info)
    {
      printf("  ERRO: %s nao tem hist_sigma.\n", nome);
      printf("  Foi gerado por uma versao antiga da macro?\n");
      return;
    }

    sigma    [i] = info->GetBinContent(1);
    sigma_err[i] = info->GetBinContent(2);
    n_aceitos[i] = info->GetBinContent(3);
    peso     [i] = info->GetBinContent(5);   // sigma / n_aceitos

    hist_sigmas->SetBinContent(i + 1, sigma[i]);
    hist_sigmas->SetBinError  (i + 1, sigma_err[i]);
    hist_sigmas->GetXaxis()->SetBinLabel(i + 1, sufixo[i]);
  }

  // ---- fator de normalizacao global
  // A maior secao de choque e' a da rodada de pThat mais baixo, mas procuro
  // em vez de supor, para a macro nao depender da ordem do vetor.
  Int_t i_maior = 0;
  for (Int_t i = 1; i < NRODADAS; i++)
    if (sigma[i] > sigma[i_maior]) i_maior = i;

  Double_t divisor = 1.;
  if (normalizar_pela_maior) divisor = peso[i_maior];

  // ---- tabela dos pesos
  printf("\n============== RODADAS E PESOS ==============\n");
  printf("  %-9s %12s %9s %9s %14s %10s\n",
         "rodada", "sigma[mb]", "err", "N", "w=sigma/N", "w/w_maior");
  for (Int_t i = 0; i < NRODADAS; i++)
    printf("  %-9s %12.4f %9.4f %9.0f %14.4e %10.5f\n",
           sufixo[i], sigma[i], sigma_err[i], n_aceitos[i],
           peso[i], peso[i] / peso[i_maior]);

  printf("---------------------------------------------\n");
  if (normalizar_pela_maior)
    printf("  normalizado pela rodada %s (a de maior sigma):\n"
           "  todo peso dividido por %.4e, entao ela fica com peso 1\n",
           sufixo[i_maior], divisor);
  else
    printf("  sem normalizacao global: eixo vertical em mb por bin\n");

  // ---- clona e aplica o peso
  for (Int_t i = 0; i < NRODADAS; i++)
  {
    TH1F *h_jato  = (TH1F*) entrada[i]->Get("hist_pt_jato");
    TH1F *h_foton = (TH1F*) entrada[i]->Get("hist_pt_foton");

    sprintf(nome, "esc_pt_jato_%s", sufixo[i]);
    esc_pt_jato[i] = (TH1F*) h_jato->Clone(nome);
    esc_pt_jato[i]->SetDirectory(arquivo);

    sprintf(nome, "esc_pt_foton_%s", sufixo[i]);
    esc_pt_foton[i] = (TH1F*) h_foton->Clone(nome);
    esc_pt_foton[i]->SetDirectory(arquivo);

    // o Scale propaga as incertezas porque o Sumw2 estava ligado na geracao
    esc_pt_jato [i]->Scale(peso[i] / divisor);
    esc_pt_foton[i]->Scale(peso[i] / divisor);

    esc_pt_jato [i]->SetLineColor(cor[i]);
    esc_pt_foton[i]->SetLineColor(cor[i]);
    esc_pt_jato [i]->SetLineWidth(2);
    esc_pt_foton[i]->SetLineWidth(2);
  }

  // ---- a soma ponderada: clona a primeira e vai somando as outras
  soma_pt_jato  = (TH1F*) esc_pt_jato [0]->Clone("soma_pt_jato");
  soma_pt_foton = (TH1F*) esc_pt_foton[0]->Clone("soma_pt_foton");
  soma_pt_jato ->SetDirectory(arquivo);
  soma_pt_foton->SetDirectory(arquivo);

  for (Int_t i = 1; i < NRODADAS; i++)
  {
    soma_pt_jato ->Add(esc_pt_jato [i]);
    soma_pt_foton->Add(esc_pt_foton[i]);
  }

  soma_pt_jato ->SetLineColor(kMagenta + 2);
  soma_pt_foton->SetLineColor(kMagenta + 2);
  soma_pt_jato ->SetLineWidth(2);
  soma_pt_foton->SetLineWidth(2);
  soma_pt_jato ->SetTitle("SOMA ponderada das 4 rodadas;pT^{jato} [GeV/c];peso por bin");
  soma_pt_foton->SetTitle("SOMA ponderada das 4 rodadas;pT^{#gamma} [GeV/c];peso por bin");

  // ---- cs1: as 4 sobrepostas, o teste de consistencia
  cs1 = new TCanvas("cs1", "4 rodadas normalizadas e sobrepostas", 1400, 600);
  cs1->Divide(2, 1);

  for (Int_t lado = 0; lado < 2; lado++)
  {
    cs1->cd(lado + 1);
    gPad->SetLogy();

    TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.88);

    for (Int_t i = 0; i < NRODADAS; i++)
    {
      TH1F *h = (lado == 0) ? esc_pt_jato[i] : esc_pt_foton[i];

      if (lado == 0)
        h->SetTitle("4 rodadas normalizadas;pT^{jato} [GeV/c];peso por bin");
      else
        h->SetTitle("4 rodadas normalizadas;pT^{#gamma} [GeV/c];peso por bin");

      if (i == 0) h->Draw("E");
      else        h->Draw("E SAME");

      sprintf(rotulo, "pThat > %s GeV", sufixo[i] + 5);
      leg->AddEntry(h, rotulo, "l");
    }
    leg->Draw();
  }
  cs1->Update();
  cs1->SaveAs("v2_combinado_sobreposto.png");

  // ---- cs2: a soma
  cs2 = new TCanvas("cs2", "Soma ponderada das 4 rodadas", 1400, 600);
  cs2->Divide(2, 1);
  cs2->cd(1); gPad->SetLogy(); soma_pt_jato ->Draw("E");
  cs2->cd(2); gPad->SetLogy(); soma_pt_foton->Draw("E");
  cs2->Update();
  cs2->SaveAs("v2_combinado_soma.png");

  // ---- cs3: as secoes de choque
  cs3 = new TCanvas("cs3", "Secoes de choque das rodadas", 800, 600);
  gPad->SetLogy();
  hist_sigmas->SetMarkerStyle(20);
  hist_sigmas->Draw("P E");
  cs3->Update();
  cs3->SaveAs("v2_combinado_sigmas.png");

  // ================================================================
  //  O TESTE: quanto cada rodada da' numa janela de pT alto, onde
  //  TODAS as 4 sao validas. Se as amostras estao certas e encaixadas,
  //  os 4 valores batem -- e a soma da' cerca de 4 vezes mais.
  // ================================================================
  Double_t pt_corte = 20.;
  Int_t    b1 = esc_pt_jato[0]->FindBin(pt_corte);
  Int_t    b2 = esc_pt_jato[0]->GetNbinsX();

  printf("=============================================\n");
  printf("  TESTE DE CONSISTENCIA\n");
  printf("  integral de pT^{jato} acima de %.0f GeV/c,\n", pt_corte);
  printf("  regiao em que as 4 rodadas sao todas validas:\n");

  Double_t soma_das_partes = 0;
  for (Int_t i = 0; i < NRODADAS; i++)
  {
    Double_t v = esc_pt_jato[i]->Integral(b1, b2);
    printf("    %-9s %12.5e\n", sufixo[i], v);
    soma_das_partes += v;
  }

  printf("  ------------------------------------------\n");
  printf("    soma das 4   %12.5e   <-- conta 4 vezes a mesma regiao\n",
         soma_das_partes);
  printf("    media das 4  %12.5e   <-- ordem de grandeza do valor fisico\n",
         soma_das_partes / NRODADAS);
  printf("  ATENCAO -- qual das 4 vale: NAO e' a de pThat mais alto.\n");
  printf("  Uma rodada de pTHatMin = X e' INCOMPLETA para observavel perto\n");
  printf("  de X: faltam os eventos de pThat abaixo de X que flutuaram para\n");
  printf("  cima. A pthat15 tem muita estatistica em pT^jato > 20, mas esta'\n");
  printf("  enviesada para BAIXO justamente por isso -- repare que ela e' a\n");
  printf("  menor das quatro.\n");
  printf("  Vale a rodada de MENOR pThat que ainda tenha estatistica aqui.\n");
  printf("  A media simples acima e' so' indicativa, nao use como resultado.\n");
  printf("\n  Se os 4 valores acima batem entre si, cada rodada esta'\n");
  printf("  normalizada certo, e a SOMA conta a mesma regiao 4 vezes.\n");
  printf("=============================================\n");

  // sem Close(): fechar destruiria os histogramas e os canvases ficariam
  // em branco. O Write() sozinho ja' grava.
  arquivo->Write();
  printf("  Arquivo salvo: jatos_foton_combinado.root\n");
}
