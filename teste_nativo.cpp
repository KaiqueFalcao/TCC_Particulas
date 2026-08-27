// Versao do teste_Ajustado.cpp usando a API NATIVA do Pythia8 em vez de TPythia8.
// Motivo: TPythia8 exige ROOT compilado com -Dpythia8=ON. A API nativa nao exige,
// e os histogramas/graficos do ROOT continuam funcionando igual.
// A fisica, os cortes e os histogramas sao os mesmos do original.

#include "Pythia8/Pythia.h"

#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"

#include <cstdio>
#include <cstdlib>

Double_t numeroBarionico(Int_t pdg) {
  Int_t a = TMath::Abs(pdg);
  if (a >= 1000 && a < 10000) return (pdg > 0) ? 1.0 : -1.0;
  return 0.0;
}

// Numero leptonico: leptons sao os codigos 11..16 (e, nu_e, mu, nu_mu, tau, nu_tau).
Double_t numeroLeptonico(Int_t pdg) {
  Int_t a = TMath::Abs(pdg);
  if (a == 11 || a == 12 || a == 13 || a == 14 || a == 15 || a == 16)
    return (pdg > 0) ? 1.0 : -1.0;
  return 0.0;
}

void distribuicoes_de_energia_e_momento(Int_t nev = 10000, Int_t ndeb = 1)
{
  Double_t soma_de_energia = 0;
  Double_t soma_px = 0, soma_py = 0, soma_pz = 0;
  Double_t px, py, pz;
  Double_t energia_part_final;
  Double_t pt, cont_pt = 0, cont_pt2 = 0;
  Double_t p_modulo, soma_modulo_p_modulo = 0;
  Double_t rapidez_part_final;

  Double_t soma_carga = 0, soma_barionico = 0, soma_leptonico = 0;
  Double_t tot_Q = 0, tot_B = 0, tot_L = 0;

  TH1F *hist_soma_energia_final = new TH1F("hist_soma_energia_final", "Soma das energias das particulas de estado final do evento", 100, 1, 1);
  TH1F *hist_momentos_particulas_finais = new TH1F("hist_pt_part_final", "Distribuicao dos momentos transversais das particulas de estado final", 100, 1, 1);
  TH1F *hist_px_part_final = new TH1F("hist_px_part_final", "Momentos transversais das particulas finais no eixo X", 100, 1, 1);
  TH1F *hist_py_part_final = new TH1F("hist_py_part_final", "Momentos transversais das particulas finais no eixo Y", 100, 1, 1);
  TH1F *hist_pz_part_final = new TH1F("hist_pz_part_final", "Momentos transversais das particulas finais no eixo Z", 100, 1, 1);
  TH2F *hist_rapidez_energia = new TH2F("hist_rapidez_energia","Energia vs rapidez;Rapidez y;Energia (GeV)",100, -15, 15, 100, 0, 200);
  TH1F *hist_soma_momento = new TH1F("hist_soma_momento", "Soma dos momentos das particulas finais dos eventos", 100, 0, 0.001);
  TH2F *dispersao_momentos = new TH2F("dispersao_momentos", "Dispersao dos momentos de estado central e final", 100, 1, 1, 100, 1, 1);
  TH1F *hist_p_particula = new TH1F("hist_p_particula", "Modulo do momento por particula", 100, 1, 1);
  TH1F *hist_soma_modulo_p = new TH1F("hist_soma_modulo_p", "Soma dos modulos por evento", 100, 1, 1);

  // --- Pythia nativo: equivalente a TPythia8 + ReadString + Initialize ---
  Pythia8::Pythia pythia;
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 43");
  pythia.readString("Beams:idA = 2212");   // Proton
  pythia.readString("Beams:idB = 2212");   // Proton
  pythia.readString("Beams:eCM = 14000.");  // 14 TeV
  if (!pythia.init()) {
    printf("ERRO: falha ao inicializar o Pythia\n");
    return;
  }

  for (Int_t iev = 0; iev < nev; iev++)
  {
    soma_modulo_p_modulo = 0;
    cont_pt = 0;
    cont_pt2 = 0;
    soma_de_energia = 0;
    soma_px = 0; soma_py = 0; soma_pz = 0;
    soma_carga = 0; soma_barionico = 0; soma_leptonico = 0;

    if (!pythia.next()) continue;

    if (iev < ndeb) pythia.event.list();

    Int_t np = pythia.event.size();

    for (Int_t ip = 0; ip < np; ip++)
    {
      const Pythia8::Particle &part = pythia.event[ip];

      Double_t eta = part.eta();

      if (part.isFinal())   // equivalente ao status > 0 do original
      {
        pt = part.pT();
        px = part.px();
        py = part.py();
        pz = part.pz();

        energia_part_final = part.e();
        rapidez_part_final = part.y();

        Int_t pdg = part.id();

        soma_px += px;
        soma_py += py;
        soma_pz += pz;

        p_modulo = TMath::Sqrt(px*px + py*py + pz*pz);
        soma_modulo_p_modulo += p_modulo;

        soma_de_energia += energia_part_final;

        // No Pythia nativo charge() ja vem em unidades de e (nao precisa dividir por 3).
        soma_carga += part.charge();
        soma_barionico += numeroBarionico(pdg);
        soma_leptonico += numeroLeptonico(pdg);

        cont_pt += 1;
        hist_momentos_particulas_finais->Fill(pt);
        hist_px_part_final->Fill(px);
        hist_py_part_final->Fill(py);
        hist_pz_part_final->Fill(pz);
        hist_p_particula->Fill(p_modulo);

        hist_rapidez_energia->Fill(rapidez_part_final, energia_part_final);

        if (-0.9 < eta && eta < 0.9) cont_pt2 += 1;
      }
    }

    dispersao_momentos->Fill(cont_pt, cont_pt2);

    Double_t modulo_p_total = TMath::Sqrt(soma_px*soma_px + soma_py*soma_py + soma_pz*soma_pz);
    hist_soma_momento->Fill(modulo_p_total);
    hist_soma_energia_final->Fill(soma_de_energia);
    hist_soma_modulo_p->Fill(soma_modulo_p_modulo);

    tot_Q += soma_carga; tot_B += soma_barionico; tot_L += soma_leptonico;
    if (iev < 5)
      printf("Evento %3d:  E=%.1f  px=%+.3f py=%+.3f pz=%+.3f  Q=%+.0f  B=%+.0f  L=%+.0f\n",
             iev, soma_de_energia, soma_px, soma_py, soma_pz,
             soma_carga, soma_barionico, soma_leptonico);
  }

  // ===== CANVAS 1: momentos por particula =====
  TCanvas *c1 = new TCanvas("c1", "Momentos por particula", 1600, 1000);
  c1->Divide(2, 2);

  c1->cd(1);
  hist_momentos_particulas_finais->SetTitle("Momento transversal (pT) por particula");
  hist_momentos_particulas_finais->GetXaxis()->SetTitle("pT (GeV)");
  hist_momentos_particulas_finais->GetYaxis()->SetTitle("Frequencia");
  hist_momentos_particulas_finais->Draw();

  c1->cd(2);
  hist_px_part_final->SetTitle("px por particula");
  hist_px_part_final->GetXaxis()->SetTitle("px (GeV)");
  hist_px_part_final->GetYaxis()->SetTitle("Frequencia");
  hist_px_part_final->Draw();

  c1->cd(3);
  hist_py_part_final->SetTitle("py por particula");
  hist_py_part_final->GetXaxis()->SetTitle("py (GeV)");
  hist_py_part_final->GetYaxis()->SetTitle("Frequencia");
  hist_py_part_final->Draw();

  c1->cd(4);
  hist_pz_part_final->SetTitle("pz por particula");
  hist_pz_part_final->GetXaxis()->SetTitle("pz (GeV)");
  hist_pz_part_final->GetYaxis()->SetTitle("Frequencia");
  hist_pz_part_final->Draw();

  c1->SaveAs("momentos_por_particula.png");

  // ===== CANVAS 2: conservacao e somas por evento =====
  TCanvas *c2 = new TCanvas("c2", "Conservacao por evento", 1600, 1000);
  c2->Divide(2, 3);

  c2->cd(1);
  hist_soma_energia_final->SetTitle("Soma da energia por evento");
  hist_soma_energia_final->GetXaxis()->SetTitle("Energia total (GeV)");
  hist_soma_energia_final->GetYaxis()->SetTitle("Frequencia");
  hist_soma_energia_final->Draw();

  c2->cd(2);
  hist_soma_momento->SetTitle("Modulo do momento total por evento");
  hist_soma_momento->GetXaxis()->SetTitle("|soma vetorial de p| (GeV)");
  hist_soma_momento->GetYaxis()->SetTitle("Frequencia");
  hist_soma_momento->Draw();

  c2->cd(3);
  hist_soma_modulo_p->SetTitle("Soma dos modulos dos momentos por evento");
  hist_soma_modulo_p->GetXaxis()->SetTitle("Soma de |p| (GeV)");
  hist_soma_modulo_p->GetYaxis()->SetTitle("Frequencia");
  hist_soma_modulo_p->Draw();

  c2->cd(4);
  dispersao_momentos->SetTitle("Multiplicidade: total x central (|eta|<0.9)");
  dispersao_momentos->GetXaxis()->SetTitle("N particulas (total)");
  dispersao_momentos->GetYaxis()->SetTitle("N particulas (central)");
  dispersao_momentos->Draw("COLZ");

  c2->cd(5);
  hist_p_particula->SetTitle("Modulo do momento por particula");
  hist_p_particula->GetXaxis()->SetTitle("|p| (GeV)");
  hist_p_particula->GetYaxis()->SetTitle("Frequencia");
  hist_p_particula->Draw();

  c2->cd(6);
  hist_rapidez_energia->SetTitle("Energia vs Rapidez");
  hist_rapidez_energia->GetXaxis()->SetTitle("Rapidez y");
  hist_rapidez_energia->GetYaxis()->SetTitle("Energia (GeV)");
  hist_rapidez_energia->Draw("COLZ");

  c2->SaveAs("conservacao_por_evento.png");

  printf("Medias sobre %d eventos:\n", nev);
  printf("  E  = %.1f GeV\n", hist_soma_energia_final->GetMean());
  printf("  px = %+.4f   py = %+.4f   pz = %+.4f GeV   (por evento)\n",
         hist_px_part_final->GetMean(), hist_py_part_final->GetMean(), hist_pz_part_final->GetMean());
  printf("  Q  = %+.3f   B  = %+.3f   L  = %+.3f\n",
         tot_Q / nev, tot_B / nev, tot_L / nev);
}

int main(int argc, char **argv)
{
  gROOT->SetBatch(kTRUE);   // nao abre janela; salva os PNGs direto
  Int_t nev  = (argc > 1) ? atoi(argv[1]) : 10000;
  Int_t ndeb = (argc > 2) ? atoi(argv[2]) : 1;
  distribuicoes_de_energia_e_momento(nev, ndeb);
  return 0;
}
