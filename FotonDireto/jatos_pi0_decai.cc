// Jatos e funcao de fragmentacao do pion neutro -- VERSAO COM O pi0 DECAINDO.
// Pythia gera, FastJet agrupa em jatos, e medimos z = pT(pi0)/pT(jato).
//
// CONJUNTOS DE HISTOGRAMAS:
//   A) z do pi0, fatiado por pT do JATO            (6 faixas)
//   B) z do pi0, fatiado por pT do PI0             (6 faixas, mesmos limites)
//   C) lider x nao-lider, sem corte                (6 faixas)
//   D) lider x nao-lider, z > corte_z              (6 faixas)
//   E) lider x nao-lider, z > corte_z e A > corte_A (6 faixas)
//
// As faixas sao as mesmas em todos: 5-8, 8-10, 10-12, 12-16, 16-20 e >20 GeV.
// A ultima e' implementada como 20 a 1000, ou seja, "acima de 20".
//
// MEDIDO com 50 mil eventos e pTHatMin = 4 (todos os pi0, sem faixa):
//   z > 0,5 -> 371 | z > 0,6 -> 175 | z > 0,7 -> 85 | z > 0,9 -> 11
//   A assimetria e' PLANA entre 0 e 1, entao A > 0,9 mantem ~10% dos pi0.

#include "TSystem.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <vector>
#include <set>

const Int_t NFAIXAS = 6;

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina, entao da' para mexer nos histogramas.
TH1F   *hist_pt_jato     = 0;
TH1F   *hist_n_part      = 0;
TH1F   *hist_pt_pi0      = 0;
TH1F   *hist_z           = 0;   // z de todos os pi0, sem fatiar
TH2F   *hist_eta_phi     = 0;
TH1F   *hist_assimetria  = 0;
TH1F   *hist_z_gama      = 0;
TH2F   *hist_z_vs_ptjato = 0;

// (A) z fatiado por pT do jato ; (B) z fatiado por pT do pi0
TH1F   *hist_z_jato[NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH1F   *hist_z_pi0 [NFAIXAS] = {0, 0, 0, 0, 0, 0};

// (C) sem corte ; (D) z > corte_z ; (E) z e A
TH2F   *hist_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F   *hist_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F   *hist_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};

TCanvas *c1 = 0;
TCanvas *c2 = 0;
TCanvas *c3 = 0;
TCanvas *c4 = 0;
TCanvas *c5 = 0;
TCanvas *c6 = 0;
TCanvas *c7 = 0;
TCanvas *c8 = 0;

Bool_t eh_neutrino(Int_t pdg) // devolve true se o PDG for de neutrino
{
  Int_t a = TMath::Abs(pdg);
  return (a == 12 || a == 14 || a == 16 || a == 18);
}

// Sobe a cadeia de maes a partir da particula i e devolve o indice do primeiro
// pi0 ancestral (-1 se nao houver).
//
// Convencao dos indices, ja' deslocada pelo ROOT no ImportParticles:
//   m1 <  0     -> nao tem mae, chegou no topo do registro
//   m2 == m1    -> "carbon copy" (mesma particula recopiada); continua subindo
//   m2 >  m1    -> FAIXA de maes = hadronizacao da string de cor. PARA AQUI.
Int_t mae_pi0(TClonesArray *lista, Int_t i)
{
  Int_t n = lista->GetEntriesFast();
  Int_t idx = i;

  for (Int_t passos = 0; passos < 200; passos++)   // trava contra laco infinito
  {
    if (idx < 0 || idx >= n) return -1;

    TParticle *p = (TParticle*) lista->At(idx);
    Int_t m1 = p->GetFirstMother();
    Int_t m2 = p->GetSecondMother();

    if (m1 < 0 || m1 >= n) return -1;   // topo do registro
    if (m2 > m1)           return -1;   // faixa de maes -> string de cor

    if (((TParticle*) lista->At(m1))->GetPdgCode() == 111) return m1;  // achou

    idx = m1;
  }
  return -1;
}

// Devolve em qual faixa um valor de pT cai (-1 se estiver fora de todas)
Int_t acha_faixa(Double_t pt, Double_t *lo, Double_t *hi)
{
  for (Int_t f = 0; f < NFAIXAS; f++)
    if (pt >= lo[f] && pt < hi[f]) return f;
  return -1;
}

void jatos_pi0_decai()
{
  // ---- parametros
  Int_t    nev          = 50000;
  Double_t R            = 0.4;
  Double_t pt_jato_min  = 5.;      // casa com a primeira faixa
  Double_t eta_part_max = 0.9;
  Double_t y_jato_max   = 0.9 - R;
  Double_t energia      = 14000.;
  Double_t corte_A      = 0.9;     // assimetria minima
  Double_t corte_z      = 0.8;     // fracao minima do pT do jato

  // seis faixas; a ultima vai ate 1000, ou seja, "acima de 20 GeV"
  Double_t faixa_lo[NFAIXAS] = { 5.,  8., 10., 12., 16.,   20.};
  Double_t faixa_hi[NFAIXAS] = { 8., 10., 12., 16., 20., 1000.};

  // ---- arquivo de saida: criado ANTES dos histogramas, para que todos
  // ---- sejam associados a ele e salvos automaticamente no Write()
  TFile *arquivo = new TFile("jatos_pi0_decai.root", "RECREATE");

  // ---- Pythia
  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("PhaseSpace:pTHatMin = 4.");  // um pouco abaixo do corte
  pythia8->ReadString("Random:seed = 43");

  pythia8->Initialize(2212, 2212, energia);

  TClonesArray *particulas = new TClonesArray("TParticle", 10000);

  // ---- FastJet
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  std::vector<fastjet::PseudoJet> entradas;

  // ---- histogramas gerais
  hist_pt_jato = new TH1F("hist_pt_jato",
      "pT dos jatos;pT^{jato} [GeV/c];jatos", 100, 0, 100);
  hist_n_part  = new TH1F("hist_n_part",
      "particulas por jato;N_{const};jatos", 60, 0, 60);
  hist_pt_pi0  = new TH1F("hist_pt_pi0",
      "pT dos #pi^{0};pT^{#pi^{0}} [GeV/c];#pi^{0}", 100, 0, 40);
  hist_z       = new TH1F("hist_z",
      "funcao de fragmentacao (todos);z = pT^{#pi^{0}}/pT^{jato};#pi^{0}",
      50, 0, 1);
  hist_eta_phi = new TH2F("hist_eta_phi",
      "mapa #eta-#phi dos #pi^{0};#eta;#phi [rad]",
      60, -eta_part_max, eta_part_max, 60, -TMath::Pi(), TMath::Pi());

  hist_assimetria = new TH1F("hist_assimetria",
      "assimetria do decaimento;A = |E_{1}-E_{2}|/(E_{1}+E_{2});#pi^{0}",
      50, 0, 1);
  hist_z_gama = new TH1F("hist_z_gama",
      "fragmentacao do foton lider;z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma",
      50, 0, 1);

  // mostra em qual regiao de pT do jato os z altos aparecem
  hist_z_vs_ptjato = new TH2F("hist_z_vs_ptjato",
      "z do #pi^{0} x pT do jato;pT^{jato} [GeV/c];z = pT^{#pi^{0}}/pT^{jato}",
      50, 0, 50, 50, 0, 1);

  // ---- histogramas fatiados por faixa.
  // Nome e titulo sao montados com sprintf porque mudam de faixa para faixa.
  // A ultima faixa recebe o rotulo "> 20" em vez de "20 - 1000".
  char nome[100], titulo[250], rotulo[50];

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "pT^{jato} > %.0f GeV", faixa_lo[f]);
    else sprintf(rotulo, "%.0f < pT^{jato} < %.0f GeV", faixa_lo[f], faixa_hi[f]);

    // (A) z fatiado por pT do JATO
    sprintf(nome, "hist_z_jato_%d", f);
    sprintf(titulo, "z do #pi^{0}, %s;z = pT^{#pi^{0}}/pT^{jato};#pi^{0}", rotulo);
    hist_z_jato[f] = new TH1F(nome, titulo, 50, 0, 1);

    // (C) lider x nao-lider, sem corte
    sprintf(nome, "hist_td_%d", f);
    sprintf(titulo, "%s, todos os #pi^{0};"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]", rotulo);
    hist_td[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);

    // (D) lider x nao-lider, so z
    sprintf(nome, "hist_Z_%d", f);
    sprintf(titulo, "%s, z > %.1f;"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]",
        rotulo, corte_z);
    hist_Z[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);

    // (E) lider x nao-lider, z e assimetria
    sprintf(nome, "hist_AZ_%d", f);
    sprintf(titulo, "%s, z > %.1f e A > %.1f;"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]",
        rotulo, corte_z, corte_A);
    hist_AZ[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);

    // (B) z fatiado por pT do PI0 -- mesmos limites de faixa
    if (f == NFAIXAS - 1) sprintf(rotulo, "pT^{#pi^{0}} > %.0f GeV", faixa_lo[f]);
    else sprintf(rotulo, "%.0f < pT^{#pi^{0}} < %.0f GeV",
                 faixa_lo[f], faixa_hi[f]);

    sprintf(nome, "hist_z_pi0_%d", f);
    sprintf(titulo, "z do #pi^{0}, %s;z = pT^{#pi^{0}}/pT^{jato};#pi^{0}", rotulo);
    hist_z_pi0[f] = new TH1F(nome, titulo, 50, 0, 1);
  }

  Long64_t n_jatos = 0, n_pi0 = 0;
  Long64_t n_gama = 0;                         // pares gama gama analisados
  Long64_t n_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_zp[NFAIXAS] = {0, 0, 0, 0, 0, 0};  // entradas do conjunto (B)

  // ---- laco de eventos
  for (Int_t iev = 0; iev < nev; iev++)
  {
    if (iev % 10000 == 0) printf("  evento %d / %d\n", iev, nev);

    pythia8->GenerateEvent();
    pythia8->ImportParticles(particulas, "All");
    Int_t np = particulas->GetEntriesFast();

    // monta a lista de entrada do FastJet
    entradas.clear();
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particulas->At(ip);

      if (part->GetStatusCode() <= 0) continue;
      if (eh_neutrino(part->GetPdgCode())) continue;
      if (part->Pt() < 1e-9) continue;
      if (TMath::Abs(part->Eta()) > eta_part_max) continue;

      fastjet::PseudoJet pj(part->Px(), part->Py(), part->Pz(), part->Energy());
      pj.set_user_index(ip);     // ponte de volta para o TClonesArray
      entradas.push_back(pj);
    }
    if (entradas.size() < 2) continue;

    // agrupa: sempre com pT minimo, senao cada particula mole vira um "jato"
    fastjet::ClusterSequence cluster_seq(entradas, jet_def);
    std::vector<fastjet::PseudoJet> jatos =
        fastjet::sorted_by_pt(cluster_seq.inclusive_jets(pt_jato_min));

    // ---- laco de jatos
    Int_t n_jatos_ev = jatos.size();
    for (Int_t ij = 0; ij < n_jatos_ev; ij++)
    {
      if (TMath::Abs(jatos[ij].rap()) > y_jato_max) continue;

      Double_t pt_jato = jatos[ij].pt();
      n_jatos++;

      std::vector<fastjet::PseudoJet> constituintes = jatos[ij].constituents();
      Int_t n_const = constituintes.size();

      hist_pt_jato->Fill(pt_jato);
      hist_n_part ->Fill(n_const);

      // faixa pelo pT do JATO (vale para os conjuntos A, C, D e E)
      Int_t faixa_jato = acha_faixa(pt_jato, faixa_lo, faixa_hi);

      // ---- laco de constituintes: procura os pi0
      std::set<Int_t> indices_pi0;

      for (Int_t ic = 0; ic < n_const; ic++)
      {
        Int_t indice = constituintes[ic].user_index();
        Int_t indice_pi0 = mae_pi0(particulas, indice);

        if (indice_pi0 < 0) continue;   // esse constituinte nao veio de pi0

        indices_pi0.insert(indice_pi0);
      }

      // agora percorremos so' os pi0 unicos encontrados acima
      std::set<Int_t>::iterator it;
      for (it = indices_pi0.begin(); it != indices_pi0.end(); it++)
      {
        TParticle *part = (TParticle*) particulas->At(*it);

        Double_t pt_pi0 = part->Pt();
        Double_t z      = pt_pi0 / pt_jato;
        n_pi0++;

        // TParticle::Phi() devolve [0, 2pi); o eixo do mapa vai de -pi a +pi
        Double_t phi = part->Phi();
        if (phi > TMath::Pi()) phi -= 2 * TMath::Pi();

        hist_pt_pi0     ->Fill(pt_pi0);
        hist_z          ->Fill(z);
        hist_eta_phi    ->Fill(part->Eta(), phi);
        hist_z_vs_ptjato->Fill(pt_jato, z);

        // (A) z fatiado por pT do JATO
        if (faixa_jato >= 0) hist_z_jato[faixa_jato]->Fill(z);

        // (B) z fatiado por pT do PI0
        Int_t faixa_pi0 = acha_faixa(pt_pi0, faixa_lo, faixa_hi);
        if (faixa_pi0 >= 0) {
          hist_z_pi0[faixa_pi0]->Fill(z);
          n_zp[faixa_pi0]++;
        }

        // O pi0 encontrado pelo mae_pi0() e' justamente a copia que decaiu,
        // entao as filhas dele sao os dois fotons.
        Int_t d1 = part->GetFirstDaughter();
        Int_t d2 = part->GetLastDaughter();

        if (d1 < 0 || d2 <= d1 || d2 >= np) continue;

        TParticle *gama1 = (TParticle*) particulas->At(d1);
        TParticle *gama2 = (TParticle*) particulas->At(d2);

        if (gama1->GetPdgCode() != 22) continue;
        if (gama2->GetPdgCode() != 22) continue;

        Double_t e1 = gama1->Energy();
        Double_t e2 = gama2->Energy();

        Double_t assimetria = TMath::Abs(e1 - e2) / (e1 + e2);
        hist_assimetria->Fill(assimetria);

        // separa o foton lider (mais energetico) do nao-lider
        TParticle *gama_lider     = gama1;
        TParticle *gama_nao_lider = gama2;
        if (e2 > e1) { gama_lider = gama2; gama_nao_lider = gama1; }

        Double_t z_gama = gama_lider->Pt() / pt_jato;
        hist_z_gama->Fill(z_gama);

        n_gama++;

        // ---- conjuntos C, D e E: lider x nao-lider por faixa de pT do jato
        if (faixa_jato >= 0)
        {
          Double_t pt_lider     = gama_lider->Pt();
          Double_t pt_nao_lider = gama_nao_lider->Pt();

          // (C) sem corte
          hist_td[faixa_jato]->Fill(pt_lider, pt_nao_lider);
          n_td[faixa_jato]++;

          // (D) so z
          if (z > corte_z) {
            hist_Z[faixa_jato]->Fill(pt_lider, pt_nao_lider);
            n_Z[faixa_jato]++;
          }

          // (E) z e assimetria
          if (z > corte_z && assimetria > corte_A) {
            hist_AZ[faixa_jato]->Fill(pt_lider, pt_nao_lider);
            n_AZ[faixa_jato]++;
          }
        }
      }
    }
  }

  // ---- canvas 1: quatro paineis gerais
  c1 = new TCanvas("c1", "Jatos e fragmentacao do pi0", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato->Draw();
  c1->cd(2); hist_n_part ->Draw();
  c1->cd(3); hist_pt_pi0 ->Draw();
  c1->cd(4); hist_z->Draw();
  c1->Update();
  c1->SaveAs("jatos_pi0_decai.png");

  // ---- canvas 2: mapa eta-phi
  c2 = new TCanvas("c2", "Mapa eta-phi dos pi0", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();
  c2->SaveAs("mapa_eta_phi_decai.png");

  // ---- canvas 3: fotons do decaimento
  c3 = new TCanvas("c3", "Fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw();                  // deve sair PLANA
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw();
  c3->Update();
  c3->SaveAs("fotons_decai.png");

  // ---- canvas 4: z versus pT do jato
  c4 = new TCanvas("c4", "z x pT do jato", 900, 700);
  gPad->SetLogz();                     // escala log nas cores: realca a cauda
  hist_z_vs_ptjato->Draw("COLZ");
  c4->Update();
  c4->SaveAs("z_vs_ptjato.png");

  // ---- canvas 5: (A) z fatiado por pT do jato
  c5 = new TCanvas("c5", "z por faixa de pT do jato", 1500, 1000);
  c5->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c5->cd(f + 1); hist_z_jato[f]->Draw(); }
  c5->Update();
  c5->SaveAs("z_por_ptjato.png");

  // ---- canvas 6: (B) z fatiado por pT do pi0
  c6 = new TCanvas("c6", "z por faixa de pT do pi0", 1500, 1000);
  c6->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c6->cd(f + 1); hist_z_pi0[f]->Draw(); }
  c6->Update();
  c6->SaveAs("z_por_ptpi0.png");

  // ---- canvas 7: (C) lider x nao-lider, sem corte
  c7 = new TCanvas("c7", "Lider x nao-lider (todos)", 1500, 1000);
  c7->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c7->cd(f + 1); hist_td[f]->Draw("COLZ"); }
  c7->Update();
  c7->SaveAs("lider_todos.png");

  // ---- canvas 8: (D) e (E) nao cabem juntos; dois arquivos separados
  c8 = new TCanvas("c8", "Lider x nao-lider (z alto)", 1500, 1000);
  c8->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c8->cd(f + 1); hist_Z[f]->Draw("COLZ"); }
  c8->Update();
  c8->SaveAs("lider_z.png");

  for (Int_t f = 0; f < NFAIXAS; f++) { c8->cd(f + 1); hist_AZ[f]->Draw("COLZ"); }
  c8->Update();
  c8->SaveAs("lider_z_e_assimetria.png");

  // ---- resumo (ANTES de fechar o arquivo: o Close apaga os histogramas)
  printf("\n===================== RESUMO =====================\n");
  printf("  eventos gerados           : %d\n", nev);
  printf("  jatos reconstruidos       : %lld\n", n_jatos);
  printf("  pi0 dentro de jatos       : %lld\n", n_pi0);
  printf("  pi0 por jato              : %.3f\n",
         n_jatos ? (Double_t) n_pi0 / n_jatos : 0.);
  printf("  <pT> dos jatos            : %.2f GeV/c\n", hist_pt_jato->GetMean());
  printf("  <z> dos pi0               : %.4f\n", hist_z->GetMean());
  printf("  pares gama gama           : %lld\n", n_gama);
  printf("  <A> assimetria            : %.4f   (esperado ~0,5 se plana)\n",
         hist_assimetria->GetMean());
  printf("  RMS da assimetria         : %.4f   (esperado ~0,289 = 1/sqrt12)\n",
         hist_assimetria->GetRMS());
  printf("--------------------------------------------------\n");
  printf("  POR FAIXA DE pT DO JATO\n");
  printf("  faixa [GeV]     todos    z>%.1f   z e A    <z>\n", corte_z);
  for (Int_t f = 0; f < NFAIXAS; f++)
    printf("  %4.0f - %-5s %8lld %8lld %8lld  %.4f\n",
           faixa_lo[f], (f == NFAIXAS - 1) ? "inf" : "",
           n_td[f], n_Z[f], n_AZ[f], hist_z_jato[f]->GetMean());
  printf("--------------------------------------------------\n");
  printf("  POR FAIXA DE pT DO PI0\n");
  printf("  faixa [GeV]   entradas     <z>\n");
  for (Int_t f = 0; f < NFAIXAS; f++)
    printf("  %4.0f - %-5s %8lld   %.4f\n",
           faixa_lo[f], (f == NFAIXAS - 1) ? "inf" : "",
           n_zp[f], hist_z_pi0[f]->GetMean());
  printf("--------------------------------------------------\n");
  printf("  pi0 acima de varios cortes de z (todas as faixas):\n");
  printf("    z > 0,5 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.5), 50));
  printf("    z > 0,6 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.6), 50));
  printf("    z > 0,7 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.7), 50));
  printf("    z > 0,8 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.8), 50));
  printf("    z > 0,9 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.9), 50));
  printf("  total de pi0 no hist_z    : %8.0f\n", hist_z->GetEntries());
  printf("==================================================\n");

  // ---- salva tudo no .root (POR ULTIMO, depois de ler os histogramas)
  arquivo->Write();
  printf("  Arquivo salvo: jatos_pi0_decai.root\n");
}