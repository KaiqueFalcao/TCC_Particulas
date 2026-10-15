// Jatos e funcao de fragmentacao do pion neutro -- VERSAO COM O pi0 DECAINDO.
// Pythia gera, FastJet agrupa em jatos, e medimos z = pT(pi0)/pT(jato).
//
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

const Int_t NFAIXAS = 5;

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina, entao da' para mexer nos histogramas.
TH1F   *hist_pt_jato = 0;
TH1F   *hist_n_part  = 0;
TH1F   *hist_pt_pi0  = 0;
TH1F   *hist_z       = 0;
TH2F   *hist_eta_phi = 0;
TH1F   *hist_assimetria = 0;
TH1F   *hist_z_gama     = 0;
TH2F   *hist_z_vs_ptjato = 0;   // NOVO: onde os z altos vivem

TH2F   *hist_lead_faixa[NFAIXAS]    = {0, 0, 0, 0, 0};   // com z > corte
TH2F   *hist_lead_faixa_td[NFAIXAS] = {0, 0, 0, 0, 0};   // todos os pi0

TCanvas *c1 = 0;
TCanvas *c2 = 0;
TCanvas *c3 = 0;
TCanvas *c4 = 0;
TCanvas *c5 = 0;
TCanvas *c6 = 0;

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

void jatos_pi0_decai()
{
  // ---- parametros
  Int_t    nev          = 50000;
  Double_t R            = 0.4;
  Double_t pt_jato_min  = 5.;      // casa com a primeira faixa
  Double_t eta_part_max = 0.9;
  Double_t y_jato_max   = 0.9 - R;
  Double_t energia      = 14000.;
  Double_t corte_z      = 0.9;     // troque aqui para 0,7 ou 0,5 se precisar

  // limites das cinco faixas de pT do jato
  Double_t faixa_lo[NFAIXAS] = { 5.,  8., 10., 12., 16.};
  Double_t faixa_hi[NFAIXAS] = { 8., 10., 12., 16., 20.};

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

  // ---- histogramas
  hist_pt_jato = new TH1F("hist_pt_jato",
      "pT dos jatos;pT^{jato} [GeV/c];jatos", 100, 0, 100);
  hist_n_part  = new TH1F("hist_n_part",
      "particulas por jato;N_{const};jatos", 60, 0, 60);
  hist_pt_pi0  = new TH1F("hist_pt_pi0",
      "pT dos #pi^{0};pT^{#pi^{0}} [GeV/c];#pi^{0}", 100, 0, 40);
  hist_z       = new TH1F("hist_z",
      "funcao de fragmentacao;z = pT^{#pi^{0}}/pT^{jato};#pi^{0}", 50, 0, 1);
  hist_eta_phi = new TH2F("hist_eta_phi",
      "mapa #eta-#phi dos #pi^{0};#eta;#phi [rad]",
      60, -eta_part_max, eta_part_max, 60, -TMath::Pi(), TMath::Pi());

  hist_assimetria = new TH1F("hist_assimetria",
      "assimetria do decaimento;A = |E_{1}-E_{2}|/(E_{1}+E_{2});#pi^{0}",
      50, 0, 1);
  hist_z_gama = new TH1F("hist_z_gama",
      "fragmentacao do foton lider;z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma",
      50, 0, 1);

  // NOVO: mostra em qual regiao de pT do jato os z altos aparecem
  hist_z_vs_ptjato = new TH2F("hist_z_vs_ptjato",
      "z do #pi^{0} x pT do jato;pT^{jato} [GeV/c];z = pT^{#pi^{0}}/pT^{jato}",
      50, 0, 50, 50, 0, 1);

  // os dez histogramas 2D (cinco faixas, com e sem corte de z).
  // Nome e titulo sao montados com sprintf porque mudam de faixa para faixa.
  char nome[100], titulo[200];

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    sprintf(nome, "hist_lead_faixa_%d", f);
    sprintf(titulo,
        "pT lider x nao-lider, %.0f < pT^{jato} < %.0f GeV, z > %.1f;"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]",
        faixa_lo[f], faixa_hi[f], corte_z);
    hist_lead_faixa[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);

    sprintf(nome, "hist_lead_faixa_td_%d", f);
    sprintf(titulo,
        "pT lider x nao-lider, %.0f < pT^{jato} < %.0f GeV, todos os #pi^{0};"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]",
        faixa_lo[f], faixa_hi[f]);
    hist_lead_faixa_td[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);
  }

  Long64_t n_jatos = 0, n_pi0 = 0;
  Long64_t n_gama = 0;                      // pares gama gama analisados
  Long64_t n_faixa[NFAIXAS]    = {0, 0, 0, 0, 0};   // entradas com z > corte
  Long64_t n_faixa_td[NFAIXAS] = {0, 0, 0, 0, 0};   // entradas sem corte

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

      // descobre em qual faixa de pT este jato cai (-1 = fora de todas)
      Int_t faixa = -1;
      for (Int_t f = 0; f < NFAIXAS; f++)
        if (pt_jato >= faixa_lo[f] && pt_jato < faixa_hi[f]) faixa = f;

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
        hist_z_vs_ptjato->Fill(pt_jato, z);   // NOVO

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

        // os histogramas 2D por faixa de pT do jato
        if (faixa >= 0)
        {
          Double_t pt_lider     = gama_lider->Pt();
          Double_t pt_nao_lider = gama_nao_lider->Pt();

          // versao sem corte de z: todos os pi0 da faixa
          hist_lead_faixa_td[faixa]->Fill(pt_lider, pt_nao_lider);
          n_faixa_td[faixa]++;

          // versao com corte: so os pi0 que levam mais que corte_z do jato
          if (z > corte_z)
          {
            hist_lead_faixa[faixa]->Fill(pt_lider, pt_nao_lider);
            n_faixa[faixa]++;
          }
        }
      }
    }
  }

  // ---- canvas 1: quatro paineis
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

  // ---- canvas 4: lider x nao-lider, com corte de z
  c4 = new TCanvas("c4", "Lider x nao-lider (com corte de z)", 1500, 1000);
  c4->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    c4->cd(f + 1);
    hist_lead_faixa[f]->Draw("COLZ");
  }
  c4->Update();
  c4->SaveAs("lider_naolider_cortez.png");

  // ---- canvas 5: lider x nao-lider, todos os pi0 (para comparar)
  c5 = new TCanvas("c5", "Lider x nao-lider (todos)", 1500, 1000);
  c5->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    c5->cd(f + 1);
    hist_lead_faixa_td[f]->Draw("COLZ");
  }
  c5->Update();
  c5->SaveAs("lider_naolider_todos.png");

  // ---- canvas 6: z versus pT do jato
  c6 = new TCanvas("c6", "z x pT do jato", 900, 700);
  gPad->SetLogz();                     // escala log nas cores: realca a cauda
  hist_z_vs_ptjato->Draw("COLZ");
  c6->Update();
  c6->SaveAs("z_vs_ptjato.png");

   // ---- resumo
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
  printf("  <z> do foton lider        : %.4f\n", hist_z_gama->GetMean());
  printf("--------------------------------------------------\n");
  printf("  ENTRADAS POR FAIXA DE pT DO JATO\n");
  printf("  faixa [GeV]     todos      z > %.1f\n", corte_z);
  for (Int_t f = 0; f < NFAIXAS; f++)
    printf("  %4.0f - %-4.0f   %8lld   %8lld\n",
           faixa_lo[f], faixa_hi[f], n_faixa_td[f], n_faixa[f]);
  printf("--------------------------------------------------\n");
  printf("  pi0 acima de varios cortes de z (todas as faixas):\n");
  printf("    z > 0,5 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.5), 50));
  printf("    z > 0,6 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.6), 50));
  printf("    z > 0,7 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.7), 50));
  printf("    z > 0,9 : %8.0f\n", hist_z->Integral(hist_z->FindBin(0.9), 50));
  printf("  total de pi0 no hist_z    : %8.0f\n", hist_z->GetEntries());
  printf("==================================================\n");

  // ---- salva tudo no arquivo .root (POR ULTIMO, depois de usar os histogramas)
  arquivo->Write();
  arquivo->Close();
  printf("  Arquivo salvo: jatos_pi0_decai.root\n");
}