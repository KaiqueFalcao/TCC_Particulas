//  ETAPA 1 -- JATOS E FUNCAO DE FRAGMENTACAO DO PION NEUTRO

//  COMPILAR
//    g++ jatos_pi0.cc -o jatos_pi0 \
//        -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8 \
//        $(fastjet-config --cxxflags --libs) $(root-config --cflags --libs)
//  RODAR
//    ./jatos_pi0
// --------------------------------------------------------------------------------

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

using namespace Pythia8;

int main()
{
  // ---- parametros 
  Int_t    nev          = 5000;    // numero de eventos
  Double_t R            = 0.4;     // raio do jato
  Double_t pt_jato_min  = 20.0;    // pT minimo do jato (GeV)
  Double_t eta_part_max = 0.9;               // aceitacao do ALICE
  Double_t y_jato_max   = 0.9 - R;           // = 0.5, para o jato caber inteiro

  // ---- configuracao do Pythia 
  Pythia pythia;
  pythia.readString("Beams:idA = 2212");            // proton
  pythia.readString("Beams:idB = 2212");            // proton
  pythia.readString("Beams:eCM = 14000.");          // 14 TeV
  pythia.readString("HardQCD:all = on");            // processos que dao jatos
  pythia.readString("PhaseSpace:pTHatMin = 20.");   // so colisoes duras
  pythia.readString("111:mayDecay = off");          // 111 = pi0, sem decair, como particula final
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 43");
  pythia.init();

  // ---- definicao do jato (FastJet) 
  fastjet::JetDefinition definicao_jato(fastjet::antikt_algorithm, R);

  // ---- histogramas 
  TH1F *hist_pt_jato = new TH1F("hist_pt_jato",
      "pT dos jatos;pT do jato (GeV);jatos", 100, 0, 300);

  TH1F *hist_n_const = new TH1F("hist_n_const",
      "Numero de particulas por jato;N particulas;jatos", 60, 0, 60);

  TH1F *hist_pt_pi0 = new TH1F("hist_pt_pi0",
      "pT dos pi0 nos jatos;pT do pi0 (GeV);pi0", 100, 0, 100);

  TH1F *hist_z = new TH1F("hist_z",
      "Funcao de fragmentacao do pi0;z_{T} = pT(pi0)/pT(jato);pi0", 50, 0, 1);

  // Mapa eta-phi dos pi0: e' o plano onde, na proxima etapa, vamos medir a
  // vizinhanca dos fotons com o raio R.
  TH2F *hist_eta_phi = new TH2F("hist_eta_phi",
      "Posicao dos pi0 no plano;eta;phi (rad)", 100, -3, 3, 100, -TMath::Pi(), TMath::Pi());

  Long64_t n_jatos = 0, n_pi0 = 0;

  //  Loop de eventos
  for (Int_t iev = 0; iev < nev; iev++)
  {
    if (!pythia.next()) continue;

    // ---- 1) monta a lista de entrada do FastJet 
    std::vector<fastjet::PseudoJet> entradas;
    // Loop de particulas, percorre todas as particulas do evento
    for (Int_t i = 0; i < pythia.event.size(); i++)
    {
      // Filtro de estado final: so particulas que sobraram no fim do evento.
      if (!pythia.event[i].isFinal()) continue;

      // Tira os neutrinos
      if (!pythia.event[i].isVisible()) continue;

      // Aceitacao: descarta o que sai colado demais na direcao do feixe.
      if (TMath::Abs(pythia.event[i].eta()) > eta_part_max) continue;

      // Quadrivetor do fastjet
      fastjet::PseudoJet p(pythia.event[i].px(), pythia.event[i].py(),
                           pythia.event[i].pz(), pythia.event[i].e());

      // Guarda o indice original: e' a ponte de volta para o registro do
      // Pythia, para depois saber QUAL particula e' cada constituinte.
      p.set_user_index(i);

      entradas.push_back(p);
    }

    if (entradas.size() < 2) continue;

    // ---- 2) agrupa em jatos 
    fastjet::ClusterSequence cs(entradas, definicao_jato);

    // inclusive_jets(pt_min) ja devolve so os jatos acima do pT minimo.
    std::vector<fastjet::PseudoJet> jatos =
        fastjet::sorted_by_pt(cs.inclusive_jets(pt_jato_min));

    // ---- 3) analisa cada jato 
    for (UInt_t j = 0; j < jatos.size(); j++)
    {
      // Aceitacao do jato em rapidez.
      if (TMath::Abs(jatos[j].rap()) > y_jato_max) continue;

      Double_t pt_jato = jatos[j].pt();
      n_jatos++;

      hist_pt_jato->Fill(pt_jato);

      // constituents() devolve as particulas que formaram este jato.
      std::vector<fastjet::PseudoJet> constituintes = jatos[j].constituents();
      hist_n_const->Fill(constituintes.size());

      // ---- 4) procura os pi0 entre os constituintes 
      for (UInt_t k = 0; k < constituintes.size(); k++)
      {
        // Volta ao registro do Pythia para saber quem e' esta particula.
        Int_t idx = constituintes[k].user_index();

        if (pythia.event[idx].id() != 111) continue;   // 111 = pi0

        Double_t pt_pi0 = constituintes[k].pt();
        Double_t z = pt_pi0 / pt_jato;                 // funcao de fragmentacao

        hist_pt_pi0->Fill(pt_pi0);
        hist_z->Fill(z);
        hist_eta_phi->Fill(pythia.event[idx].eta(), pythia.event[idx].phi());

        n_pi0++;
      }
    }

    if (iev < 5)
      printf("Evento %2d:  %2d jatos reconstruidos\n", iev, (Int_t)jatos.size());
  }

  //  RESUMO
  printf("\n""RESULTADOS:\n");
  printf("  eventos gerados     : %d\n", nev);
  printf("  jatos reconstruidos : %lld\n", n_jatos);
  printf("  pi0 dentro de jatos : %lld\n", n_pi0);
  printf("  pi0 por jato        : %.2f\n",
         n_jatos > 0 ? (Double_t)n_pi0 / n_jatos : 0.0);
  printf("  <pT> dos jatos      : %.1f GeV\n", hist_pt_jato->GetMean());
  printf("  <z_T> dos pi0       : %.3f\n", hist_z->GetMean());

  //  GRAFICOS
  TCanvas *c1 = new TCanvas("c1", "Jatos e fragmentacao", 1400, 1000);
  c1->Divide(2, 2);

  c1->cd(1);
  hist_pt_jato->Draw();

  c1->cd(2);
  hist_n_const->Draw();

  c1->cd(3);
  hist_pt_pi0->Draw();

  c1->cd(4);
  hist_z->Draw();

  c1->SaveAs("jatos_pi0.png");

  TCanvas *c2 = new TCanvas("c2", "Mapa eta-phi dos pi0", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->SaveAs("mapa_eta_phi.png");

  return 0;
}