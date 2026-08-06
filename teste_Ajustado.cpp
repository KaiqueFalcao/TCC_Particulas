
#include "TSystem.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h" // ADICIONADO
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"

// ADICIONADO
Double_t numeroBarionico(Int_t pdg) {
  Int_t a = TMath::Abs(pdg);
  if (a >= 1000 && a < 10000) return (pdg > 0) ? 1.0 : -1.0;
  return 0.0;
}
// Numero leptonico: leptons sao os codigos 11..16 (e, nu_e, mu, nu_mu, tau, nu_tau).
Double_t numeroLeptonico(Int_t pdg) {
  Int_t a = TMath::Abs(pdg);
  if (a == 11 /*eletron */|| a == 12 /*nu_eletron */|| a == 13 /*muon */ || a == 14 /*nu_mon*/ || a == 15 /*tau */ || a == 16 /*nu_tau */)
    return (pdg > 0) ? 1.0 : -1.0;
  return 0.0;
}

void distribuicoes_de_energia_e_momento(Int_t nev = 10000, Int_t ndeb = 1 /* Listagem */)
{
  gSystem->Load("libEG");
  gSystem->Load("libEgPythia8");


  Double_t soma_de_energia=0; // CORRIGIDO
  Double_t soma_px = 0, soma_py = 0, soma_pz = 0;     // CORRIGIDO (soma vetorial)
  Double_t px, py, pz; // ADICIONADO (Faltou pz)
  Double_t energia_part_final;
  Double_t pt, cont_pt = 0, cont_pt2 = 0;
  Double_t p_modulo, soma_modulo_p_modulo = 0;
  Double_t rapidez_part_final;
  

  // ADICIONADO: acumuladores dos numeros quanticos que faltavam no original.
  Double_t soma_carga = 0, soma_barionico = 0, soma_leptonico = 0;
  Double_t tot_Q = 0, tot_B = 0, tot_L = 0;   // para tirar a media no fim

  TH1F *hist_soma_energia_final = new TH1F("hist_soma_energia_final", "Soma das energias das particulas de estado final do evento", 100, 1, 1);
  TH1F *hist_momentos_particulas_finais = new TH1F("hist_pt_part_final", "Distribuicao dos momentos transversais das particulas de estado final", 100, 1, 1);
  TH1F *hist_px_part_final = new TH1F("hist_px_part_final", "Momentos transversais das particulas finais no eixo X", 100, 1, 1);
  TH1F *hist_py_part_final = new TH1F("hist_py_part_final", "Momentos transversais das particulas finais no eixo Y", 100, 1, 1);
  TH1F *hist_pz_part_final = new TH1F("hist_pz_part_final", "Momentos transversais das particulas finais no eixo Z", 100, 1, 1);
  TH2F *hist_rapidez_energia = new TH2F("hist_rapidez_energia","Energia vs rapidez;Rapidez y;Energia (GeV)",100, -15, 15, 100, 0, 200);      // eixo y: energia, de 0 a 200
  TH1F *hist_soma_momento = new TH1F("hist_soma_momento", "Soma dos momentos das particulas finais dos eventos", 100, 0, 0.001);
  TH2F *dispersao_momentos = new TH2F("dispersao_momentos", "Dispersao dos momentos de estado central e final", 100, 1, 1, 100, 1, 1);
  TH1F *hist_p_particula = new TH1F("hist_p_particula", "Modulo do momento por particula", 100, 1, 1);
  TH1F *hist_soma_modulo_p = new TH1F("hist_soma_modulo_p", "Soma dos modulos por evento", 100, 1, 1);

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 43");

  pythia8->Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* GeV */); /* 14 TeV = 14000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 10000);

  for (Int_t iev = 0; iev < nev; iev++)
  {
    soma_modulo_p_modulo = 0;
    cont_pt = 0;
    cont_pt2 = 0;
    soma_de_energia = 0;                                       // CORRIGIDO (faltava zerar a energia)
    soma_px = 0; soma_py = 0; soma_pz = 0;                     // CORRIGIDO
    soma_carga = 0; soma_barionico = 0; soma_leptonico = 0;    // ADICIONADO

    pythia8->GenerateEvent();

    if (iev < ndeb)
      pythia8->EventListing();

    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle *)particles->At(ip); /* Criacao do ponteiro */

      Int_t status = part->GetStatusCode();
      Double_t eta = part->Eta();
      Double_t phi = part->Phi();

      if (status > 0) /* Particula de estado final */
      {
        pt = part->Pt();
        px = part->Px();
        py = part->Py();
        pz = part->Pz();                     // ADICIONADO

        energia_part_final = part->Energy();
        rapidez_part_final = part->Y();

        Int_t pdg = part->GetPdgCode();      // ADICIONADO: identidade da particula

        // CORRIGIDO: momento conservado = soma vetorial das componentes.
        soma_px += px;
        soma_py += py;
        soma_pz += pz;

        p_modulo = TMath::Sqrt(px*px + py*py + pz*pz); // ADICIONADO2
        soma_modulo_p_modulo += p_modulo;             // ADICIONADO2

        soma_de_energia += energia_part_final;

        // ADICIONADO: os numeros quanticos que faltavam.
        // Carga: o ROOT devolve em unidades de e/3
        TParticlePDG *info = part->GetPDG();
        if (info) soma_carga += info->Charge() / 3.0;
        soma_barionico += numeroBarionico(pdg);
        soma_leptonico += numeroLeptonico(pdg);

        cont_pt += 1;
        hist_momentos_particulas_finais->Fill(pt);
        hist_px_part_final->Fill(px);
        hist_py_part_final->Fill(py);
        hist_pz_part_final->Fill(pz); // ADICIONADO: histograma de pz
        hist_p_particula->Fill(p_modulo);

        hist_rapidez_energia->Fill(rapidez_part_final, energia_part_final); // ADICIONADO: histograma de rapidez

        if (-0.9 < eta && eta < 0.9)
        {
          cont_pt2 += 1;
        }
      }
    }
    
    dispersao_momentos->Fill(cont_pt, cont_pt2);

    Double_t modulo_p_total = TMath::Sqrt(soma_px*soma_px + soma_py*soma_py + soma_pz*soma_pz);
    hist_soma_momento->Fill(modulo_p_total);   // CORRIGIDO
    hist_soma_energia_final->Fill(soma_de_energia); // CORRIGIDO (Posição dentro do loop)
    hist_soma_modulo_p->Fill(soma_modulo_p_modulo); // ADICIONADO: histograma da soma dos modulos por evento

    // ADICIONADO: acumula para a media e imprime os primeiros eventos.
    tot_Q += soma_carga; tot_B += soma_barionico; tot_L += soma_leptonico;
    if (iev < 5)
      printf("Evento %3d:  E=%.1f  px=%+.3f py=%+.3f pz=%+.3f  Q=%+.0f  B=%+.0f  L=%+.0f\n",
             iev, soma_de_energia, soma_px, soma_py, soma_pz,
             soma_carga, soma_barionico, soma_leptonico);

  }


 // ===== CANVAS 1: momentos por particula (retratos) =====
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

  // ADICIONADO
  printf("Medias sobre %d eventos:\n", nev);
  printf("  E  = %.1f GeV\n", hist_soma_energia_final->GetMean());
  printf("  px = %+.4f   py = %+.4f   pz = %+.4f GeV   (por evento)\n",
         hist_px_part_final->GetMean(), hist_py_part_final->GetMean(), hist_pz_part_final->GetMean());
  printf("  Q  = %+.3f   B  = %+.3f   L  = %+.3f\n",
         tot_Q / nev, tot_B / nev, tot_L / nev);


}