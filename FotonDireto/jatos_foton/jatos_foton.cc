// Jatos e funcao de fragmentacao do FOTON -- copia do jatos_pi0_decai.cc.
// Pythia gera, FastJet agrupa em jatos, e medimos z = pT(foton)/pT(jato).
//
// Diferenca para o jatos_pi0_decai.cc: o objeto estudado deixa de ser o pi0
// (achado subindo a cadeia de maes) e passa a ser o FOTON (PDG 22), que e'
// particula de estado final e entra direto como constituinte do jato.
// Para cada foton candidato marcamos:
//   - origem     : a mae direta e' um pi0 (111)?             -> veio_de_pi0
//   - isolamento : soma do pT num cone R_iso em volta dele    -> isolado
//
// CONJUNTOS DE HISTOGRAMAS:
//   A) z do foton, fatiado por pT do JATO              (6 faixas)
//   B) z do foton, fatiado por pT do FOTON             (6 faixas, mesmos limites)
//   C) lider x nao-lider dos pi0, sem corte            (6 faixas)  -- como antes
//   D) lider x nao-lider dos pi0, z > corte_z          (6 faixas)  -- como antes
//   E) lider x nao-lider dos pi0, z > corte_z e A      (6 faixas)  -- como antes
//   F) z do foton nos 4 grupos isolamento x origem     (6 faixas de pT do jato)
//
// As faixas sao as mesmas em todos: 5-8, 8-10, 10-12, 12-16, 16-20 e >20 GeV.
// A ultima e' implementada como 20 a 1000, ou seja, "acima de 20".

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
TH1F   *hist_pt_foton    = 0;
TH1F   *hist_z           = 0;   // z de todos os fotons, sem fatiar
TH2F   *hist_eta_phi     = 0;
TH1F   *hist_assimetria  = 0;
TH1F   *hist_z_gama      = 0;   // foton LIDER do decaimento do pi0 (como antes)
TH2F   *hist_z_vs_ptjato = 0;
TH1F   *hist_soma_cone   = 0;   // soma de pT no cone de isolamento

// (A) z fatiado por pT do jato ; (B) z fatiado por pT do foton
TH1F   *hist_z_jato [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH1F   *hist_z_foton[NFAIXAS] = {0, 0, 0, 0, 0, 0};

// (C) sem corte ; (D) z > corte_z ; (E) z e A
TH2F   *hist_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F   *hist_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F   *hist_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};

// (F) z do foton nos quatro grupos isolamento x origem
TH1F   *hist_z_iso_pi0       [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH1F   *hist_z_iso_nao_pi0   [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH1F   *hist_z_naoiso_pi0    [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH1F   *hist_z_naoiso_nao_pi0[NFAIXAS] = {0, 0, 0, 0, 0, 0};

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
TCanvas *c13 = 0;
TCanvas *c14 = 0;

Bool_t eh_neutrino(Int_t pdg) // devolve true se o PDG for de neutrino
{
  Int_t a = TMath::Abs(pdg);
  return (a == 12 || a == 14 || a == 16 || a == 18);
}

// Soma do pT das particulas finais num cone de raio R_iso em volta do foton
// de indice i_foton. Percorre o evento INTEIRO, nao so' o jato: o isolamento
// olha tudo que esta' em volta do foton, esteja ou nao no mesmo jato.
Double_t soma_cone(TClonesArray *lista, Int_t i_foton,
                   Double_t R_iso, Double_t eta_part_max)
{
  Int_t n = lista->GetEntriesFast();

  TParticle *foton  = (TParticle*) lista->At(i_foton);
  Double_t eta_foton = foton->Eta();
  Double_t phi_foton = foton->Phi();

  Double_t soma_pt_cone = 0;

  for (Int_t i = 0; i < n; i++)
  {
    if (i == i_foton) continue;                       // o proprio foton

    TParticle *p = (TParticle*) lista->At(i);

    if (p->GetStatusCode() <= 0) continue;            // intermediaria: contaria a mesma energia duas vezes
    if (eh_neutrino(p->GetPdgCode())) continue;
    if (TMath::Abs(p->Eta()) > eta_part_max) continue;

    Double_t d_eta = p->Eta() - eta_foton;

    // Phi() devolve [0, 2pi); a diferenca e' trazida para [-pi, pi]
    Double_t d_phi = p->Phi() - phi_foton;
    if (d_phi >  TMath::Pi()) d_phi -= 2 * TMath::Pi();
    if (d_phi < -TMath::Pi()) d_phi += 2 * TMath::Pi();

    Double_t R_i = TMath::Sqrt(d_eta * d_eta + d_phi * d_phi);
    if (R_i > R_iso) continue;

    soma_pt_cone += p->Pt();
  }
  return soma_pt_cone;
}

// Devolve em qual faixa um valor de pT cai (-1 se estiver fora de todas)
Int_t acha_faixa(Double_t pt, Double_t *lo, Double_t *hi)
{
  for (Int_t f = 0; f < NFAIXAS; f++)
    if (pt >= lo[f] && pt < hi[f]) return f;
  return -1;
}

void jatos_foton()
{
  // ---- parametros
  Int_t    nev          = 10000;   // teste
  Double_t R            = 0.4;
  Double_t pt_jato_min  = 5.;      // casa com a primeira faixa
  Double_t eta_part_max = 0.9;
  Double_t y_jato_max   = 0.9 - R;
  Double_t energia      = 14000.;
  Double_t corte_A      = 0.9;     // assimetria minima
  Double_t corte_z      = 0.8;     // fracao minima do pT do jato
  Double_t R_iso        = 0.4;     // raio do cone de isolamento
  Double_t pt_iso_max   = 2.0;     // GeV/c: isolado se a soma no cone for menor

  // seis faixas; a ultima vai ate 1000, ou seja, "acima de 20 GeV"
  Double_t faixa_lo[NFAIXAS] = { 5.,  8., 10., 12., 16.,   20.};
  Double_t faixa_hi[NFAIXAS] = { 8., 10., 12., 16., 20., 1000.};

  // ---- arquivo de saida: criado ANTES dos histogramas, para que todos
  // ---- sejam associados a ele e salvos automaticamente no Write()
  TFile *arquivo = new TFile("jatos_foton.root", "RECREATE");

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
  char nome[100], titulo[250], rotulo[50];

  hist_pt_jato  = new TH1F("hist_pt_jato",
      "pT dos jatos;pT^{jato} [GeV/c];jatos", 100, 0, 100);
  hist_n_part   = new TH1F("hist_n_part",
      "particulas por jato;N_{const};jatos", 60, 0, 60);
  hist_pt_foton = new TH1F("hist_pt_foton",
      "pT dos #gamma;pT^{#gamma} [GeV/c];#gamma", 100, 0, 40);
  hist_z        = new TH1F("hist_z",
      "funcao de fragmentacao do foton (todos);z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma",
      50, 0, 1);
  hist_eta_phi  = new TH2F("hist_eta_phi",
      "mapa #eta-#phi dos #gamma;#eta;#phi [rad]",
      60, -eta_part_max, eta_part_max, 60, -TMath::Pi(), TMath::Pi());

  hist_assimetria = new TH1F("hist_assimetria",
      "assimetria do decaimento;A = |E_{1}-E_{2}|/(E_{1}+E_{2});#pi^{0}",
      50, 0, 1);
  hist_z_gama = new TH1F("hist_z_gama",
      "fragmentacao do foton lider;z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma",
      50, 0, 1);

  // mostra em qual regiao de pT do jato os z altos aparecem
  hist_z_vs_ptjato = new TH2F("hist_z_vs_ptjato",
      "z do #gamma x pT do jato;pT^{jato} [GeV/c];z_{#gamma} = pT^{#gamma}/pT^{jato}",
      50, 0, 50, 50, 0, 1);

  // isolamento: 0,5 GeV por bin, entao o corte de 2 GeV cai numa borda de bin
  sprintf(titulo, "soma de pT no cone de isolamento (R < %.1f);"
      "#Sigma pT^{cone} [GeV/c];#gamma", R_iso);
  hist_soma_cone = new TH1F("hist_soma_cone", titulo, 80, 0, 40);

  // ---- histogramas fatiados por faixa.
  // Nome e titulo sao montados com sprintf porque mudam de faixa para faixa.
  // A ultima faixa recebe o rotulo "> 20" em vez de "20 - 1000".
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "pT^{jato} > %.0f GeV", faixa_lo[f]);
    else sprintf(rotulo, "%.0f < pT^{jato} < %.0f GeV", faixa_lo[f], faixa_hi[f]);

    // (A) z fatiado por pT do JATO
    sprintf(nome, "hist_z_jato_%d", f);
    sprintf(titulo, "z do #gamma, %s;z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", rotulo);
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

    // (F) z do foton nos quatro grupos isolamento x origem
    sprintf(nome, "hist_z_iso_pi0_%d", f);
    sprintf(titulo, "#gamma isolado e de #pi^{0}, %s;"
        "z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", rotulo);
    hist_z_iso_pi0[f] = new TH1F(nome, titulo, 50, 0, 1);

    sprintf(nome, "hist_z_iso_nao_pi0_%d", f);
    sprintf(titulo, "#gamma isolado e nao de #pi^{0}, %s;"
        "z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", rotulo);
    hist_z_iso_nao_pi0[f] = new TH1F(nome, titulo, 50, 0, 1);

    sprintf(nome, "hist_z_naoiso_pi0_%d", f);
    sprintf(titulo, "#gamma nao isolado e de #pi^{0}, %s;"
        "z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", rotulo);
    hist_z_naoiso_pi0[f] = new TH1F(nome, titulo, 50, 0, 1);

    sprintf(nome, "hist_z_naoiso_nao_pi0_%d", f);
    sprintf(titulo, "#gamma nao isolado e nao de #pi^{0}, %s;"
        "z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", rotulo);
    hist_z_naoiso_nao_pi0[f] = new TH1F(nome, titulo, 50, 0, 1);

    // (B) z fatiado por pT do FOTON -- mesmos limites de faixa
    if (f == NFAIXAS - 1) sprintf(rotulo, "pT^{#gamma} > %.0f GeV", faixa_lo[f]);
    else sprintf(rotulo, "%.0f < pT^{#gamma} < %.0f GeV",
                 faixa_lo[f], faixa_hi[f]);

    sprintf(nome, "hist_z_foton_%d", f);
    sprintf(titulo, "z do #gamma, %s;z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", rotulo);
    hist_z_foton[f] = new TH1F(nome, titulo, 50, 0, 1);
  }

  Long64_t n_jatos = 0, n_fotons = 0;
  Long64_t n_fotons_pi0 = 0, n_isolados = 0;
  Long64_t n_gama = 0;                         // pares gama gama analisados
  Long64_t n_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_zf[NFAIXAS] = {0, 0, 0, 0, 0, 0};  // entradas do conjunto (B)

  // contagens do conjunto (F), por faixa de pT do jato
  Long64_t n_iso_pi0       [NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_iso_nao_pi0   [NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_naoiso_pi0    [NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_naoiso_nao_pi0[NFAIXAS] = {0, 0, 0, 0, 0, 0};

  // ---- laco de eventos
  for (Int_t iev = 0; iev < nev; iev++)
  {
    if (iev % 2000 == 0) printf("  evento %d / %d\n", iev, nev);

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

      // faixa pelo pT do JATO (vale para os conjuntos A, C, D, E e F)
      Int_t faixa_jato = acha_faixa(pt_jato, faixa_lo, faixa_hi);

      // pi0 de onde vieram os fotons deste jato. O set so' serve para nao
      // repetir o pi0 na assimetria quando os DOIS fotons dele estao no jato;
      // os candidatos em si sao achados direto pelo PDG, sem o set.
      std::set<Int_t> indices_pi0;

      // ---- laco de constituintes: os candidatos sao os fotons
      for (Int_t ic = 0; ic < n_const; ic++)
      {
        Int_t indice = constituintes[ic].user_index();
        TParticle *foton = (TParticle*) particulas->At(indice);

        if (foton->GetPdgCode() != 22) continue;   // so' fotons

        Double_t pt_foton = foton->Pt();
        Double_t z_gama   = pt_foton / pt_jato;
        n_fotons++;

        // ---- origem: a mae direta e' um pi0?
        // GetFirstMother() devolve um INDICE na lista, nao um PDG.
        Int_t mae = foton->GetFirstMother();
        Bool_t veio_de_pi0 = kFALSE;
        if (mae >= 0 && mae < np)
          veio_de_pi0 = (((TParticle*) particulas->At(mae))->GetPdgCode() == 111);

        if (veio_de_pi0) {
          n_fotons_pi0++;
          indices_pi0.insert(mae);
        }

        // ---- isolamento
        Double_t soma_pt_cone = soma_cone(particulas, indice, R_iso, eta_part_max);
        Bool_t isolado = (soma_pt_cone < pt_iso_max);
        if (isolado) n_isolados++;

        // TParticle::Phi() devolve [0, 2pi); o eixo do mapa vai de -pi a +pi
        Double_t phi = foton->Phi();
        if (phi > TMath::Pi()) phi -= 2 * TMath::Pi();

        hist_pt_foton   ->Fill(pt_foton);
        hist_z          ->Fill(z_gama);
        hist_eta_phi    ->Fill(foton->Eta(), phi);
        hist_z_vs_ptjato->Fill(pt_jato, z_gama);
        hist_soma_cone  ->Fill(soma_pt_cone);

        if (faixa_jato >= 0)
        {
          // (A) z fatiado por pT do JATO
          hist_z_jato[faixa_jato]->Fill(z_gama);

          // (F) quatro grupos isolamento x origem
          if (isolado && veio_de_pi0) {
            hist_z_iso_pi0[faixa_jato]->Fill(z_gama);
            n_iso_pi0[faixa_jato]++;
          }
          else if (isolado && !veio_de_pi0) {
            hist_z_iso_nao_pi0[faixa_jato]->Fill(z_gama);
            n_iso_nao_pi0[faixa_jato]++;
          }
          else if (!isolado && veio_de_pi0) {
            hist_z_naoiso_pi0[faixa_jato]->Fill(z_gama);
            n_naoiso_pi0[faixa_jato]++;
          }
          else {
            hist_z_naoiso_nao_pi0[faixa_jato]->Fill(z_gama);
            n_naoiso_nao_pi0[faixa_jato]++;
          }
        }

        // (B) z fatiado por pT do FOTON
        Int_t faixa_foton = acha_faixa(pt_foton, faixa_lo, faixa_hi);
        if (faixa_foton >= 0) {
          hist_z_foton[faixa_foton]->Fill(z_gama);
          n_zf[faixa_foton]++;
        }
      }

      // ---- assimetria e lider x nao-lider: igual ao jatos_pi0_decai.cc,
      // ---- para cada pi0 (sem repetir) que deixou pelo menos um foton no jato
      std::set<Int_t>::iterator it;
      for (it = indices_pi0.begin(); it != indices_pi0.end(); it++)
      {
        TParticle *part = (TParticle*) particulas->At(*it);

        Double_t z_pi0 = part->Pt() / pt_jato;   // usado nos cortes D e E, como antes

        // A mae direta do foton e' a copia do pi0 que decaiu,
        // entao as filhas dela sao os dois fotons.
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
          if (z_pi0 > corte_z) {
            hist_Z[faixa_jato]->Fill(pt_lider, pt_nao_lider);
            n_Z[faixa_jato]++;
          }

          // (E) z e assimetria
          if (z_pi0 > corte_z && assimetria > corte_A) {
            hist_AZ[faixa_jato]->Fill(pt_lider, pt_nao_lider);
            n_AZ[faixa_jato]++;
          }
        }
      }
    }
  }

  // ---- canvas 1: quatro paineis gerais
  c1 = new TCanvas("c1", "Jatos e fragmentacao do foton", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato ->Draw();
  c1->cd(2); hist_n_part  ->Draw();
  c1->cd(3); hist_pt_foton->Draw();
  c1->cd(4); hist_z->Draw();
  c1->Update();
  c1->SaveAs("foton_jatos.png");

  // ---- canvas 2: mapa eta-phi
  c2 = new TCanvas("c2", "Mapa eta-phi dos fotons", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();
  c2->SaveAs("foton_mapa_eta_phi.png");

  // ---- canvas 3: fotons do decaimento do pi0 (como antes)
  c3 = new TCanvas("c3", "Fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw();                  // deve sair PLANA
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw();
  c3->Update();
  c3->SaveAs("foton_assimetria.png");

  // ---- canvas 4: z versus pT do jato
  c4 = new TCanvas("c4", "z x pT do jato", 900, 700);
  gPad->SetLogz();                     // escala log nas cores: realca a cauda
  hist_z_vs_ptjato->Draw("COLZ");
  c4->Update();
  c4->SaveAs("foton_z_vs_ptjato.png");

  // ---- canvas 5: (A) z fatiado por pT do jato
  c5 = new TCanvas("c5", "z do foton por faixa de pT do jato", 1500, 1000);
  c5->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c5->cd(f + 1); hist_z_jato[f]->Draw(); }
  c5->Update();
  c5->SaveAs("foton_z_por_ptjato.png");

  // ---- canvas 6: (B) z fatiado por pT do foton
  c6 = new TCanvas("c6", "z do foton por faixa de pT do foton", 1500, 1000);
  c6->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c6->cd(f + 1); hist_z_foton[f]->Draw(); }
  c6->Update();
  c6->SaveAs("foton_z_por_ptfoton.png");

  // ---- canvas 7: (C) lider x nao-lider, sem corte
  c7 = new TCanvas("c7", "Lider x nao-lider (todos)", 1500, 1000);
  c7->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c7->cd(f + 1); hist_td[f]->Draw("COLZ"); }
  c7->Update();
  c7->SaveAs("foton_lider_todos.png");

  // ---- canvas 8: (D) lider x nao-lider, z alto
  c8 = new TCanvas("c8", "Lider x nao-lider (z alto)", 1500, 1000);
  c8->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c8->cd(f + 1); hist_Z[f]->Draw("COLZ"); }
  c8->Update();
  c8->SaveAs("foton_lider_z.png");

  // ---- canvas 9: (E) lider x nao-lider, z alto e A alta
  c9 = new TCanvas("c9", "Lider x nao-lider (z e A altos)", 1500, 1000);
  c9->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c9->cd(f + 1); hist_AZ[f]->Draw("COLZ"); }
  c9->Update();
  c9->SaveAs("foton_lider_z_e_assimetria.png");

  // ---- canvas 10: soma de pT no cone de isolamento
  c10 = new TCanvas("c10", "Isolamento dos fotons", 900, 700);
  gPad->SetLogy();
  hist_soma_cone->Draw();
  c10->Update();
  c10->SaveAs("foton_soma_cone.png");

  // ---- canvas 11 a 14: (F) z do foton nos quatro grupos
  c11 = new TCanvas("c11", "z: isolado e de pi0", 1500, 1000);
  c11->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c11->cd(f + 1); hist_z_iso_pi0[f]->Draw(); }
  c11->Update();
  c11->SaveAs("foton_z_iso_pi0.png");

  c12 = new TCanvas("c12", "z: isolado e nao de pi0", 1500, 1000);
  c12->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c12->cd(f + 1); hist_z_iso_nao_pi0[f]->Draw(); }
  c12->Update();
  c12->SaveAs("foton_z_iso_nao_pi0.png");

  c13 = new TCanvas("c13", "z: nao isolado e de pi0", 1500, 1000);
  c13->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c13->cd(f + 1); hist_z_naoiso_pi0[f]->Draw(); }
  c13->Update();
  c13->SaveAs("foton_z_naoiso_pi0.png");

  c14 = new TCanvas("c14", "z: nao isolado e nao de pi0", 1500, 1000);
  c14->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c14->cd(f + 1); hist_z_naoiso_nao_pi0[f]->Draw(); }
  c14->Update();
  c14->SaveAs("foton_z_naoiso_nao_pi0.png");

  // ---- resumo (ANTES de fechar o arquivo: o Close apaga os histogramas)
  printf("\n===================== RESUMO =====================\n");
  printf("  eventos gerados           : %d\n", nev);
  printf("  jatos reconstruidos       : %lld\n", n_jatos);
  printf("  fotons dentro de jatos    : %lld\n", n_fotons);
  printf("  fotons por jato           : %.3f\n",
         n_jatos ? (Double_t) n_fotons / n_jatos : 0.);
  printf("  <pT> dos jatos            : %.2f GeV/c\n", hist_pt_jato->GetMean());
  printf("  <z> dos fotons            : %.4f\n", hist_z->GetMean());
  printf("  fotons vindos de pi0      : %lld  (%.1f%%)\n", n_fotons_pi0,
         n_fotons ? 100. * n_fotons_pi0 / n_fotons : 0.);
  printf("  fotons isolados           : %lld  (%.1f%%)\n", n_isolados,
         n_fotons ? 100. * n_isolados / n_fotons : 0.);
  printf("  pares gama gama           : %lld\n", n_gama);
  printf("  <A> assimetria            : %.4f   (esperado ~0,5 se plana)\n",
         hist_assimetria->GetMean());
  printf("  RMS da assimetria         : %.4f   (esperado ~0,289 = 1/sqrt12)\n",
         hist_assimetria->GetRMS());
  printf("--------------------------------------------------\n");
  printf("  POR FAIXA DE pT DO JATO (lider x nao-lider dos pi0)\n");
  printf("  faixa [GeV]     todos  z_pi0>%.1f  z e A   <z_gama>\n", corte_z);
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "> %.0f", faixa_lo[f]);
    else sprintf(rotulo, "%.0f - %.0f", faixa_lo[f], faixa_hi[f]);
    printf("  %-11s %9lld %9lld %8lld    %.4f\n",
           rotulo, n_td[f], n_Z[f], n_AZ[f], hist_z_jato[f]->GetMean());
  }
  printf("--------------------------------------------------\n");
  printf("  POR FAIXA DE pT DO FOTON\n");
  printf("  faixa [GeV]   entradas     <z_gama>\n");
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "> %.0f", faixa_lo[f]);
    else sprintf(rotulo, "%.0f - %.0f", faixa_lo[f], faixa_hi[f]);
    printf("  %-11s %9lld     %.4f\n", rotulo, n_zf[f], hist_z_foton[f]->GetMean());
  }
  printf("--------------------------------------------------\n");
  printf("  fotons acima de varios cortes de z (todas as faixas):\n");
  printf("    z > 0,5 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.5), 50));
  printf("    z > 0,6 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.6), 50));
  printf("    z > 0,7 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.7), 50));
  printf("    z > 0,8 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.8), 50));
  printf("    z > 0,9 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.9), 50));
  printf("  total de fotons no hist_z : %8.0f\n", hist_z->GetEntries());
  printf("--------------------------------------------------\n");
  printf("  ISOLAMENTO x ORIGEM, POR FAIXA DE pT DO JATO\n");
  printf("  (isolado: soma de pT no cone R < %.1f menor que %.1f GeV/c)\n",
         R_iso, pt_iso_max);
  printf("  faixa [GeV]    iso+pi0  iso+outro  naoiso+pi0  naoiso+outro\n");

  Long64_t t_iso_pi0 = 0, t_iso_nao_pi0 = 0, t_naoiso_pi0 = 0, t_naoiso_nao_pi0 = 0;
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "> %.0f", faixa_lo[f]);
    else sprintf(rotulo, "%.0f - %.0f", faixa_lo[f], faixa_hi[f]);
    printf("  %-11s %10lld %10lld  %10lld    %10lld\n", rotulo,
           n_iso_pi0[f], n_iso_nao_pi0[f], n_naoiso_pi0[f], n_naoiso_nao_pi0[f]);

    t_iso_pi0        += n_iso_pi0[f];
    t_iso_nao_pi0    += n_iso_nao_pi0[f];
    t_naoiso_pi0     += n_naoiso_pi0[f];
    t_naoiso_nao_pi0 += n_naoiso_nao_pi0[f];
  }
  printf("  %-11s %10lld %10lld  %10lld    %10lld\n", "total",
         t_iso_pi0, t_iso_nao_pi0, t_naoiso_pi0, t_naoiso_nao_pi0);
  printf("==================================================\n");

  // ---- salva tudo no .root (POR ULTIMO, depois de ler os histogramas)
  arquivo->Write();
  arquivo->Close();
  printf("  Arquivo salvo: jatos_foton.root\n");
}
