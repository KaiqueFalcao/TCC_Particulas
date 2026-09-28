// ===========================================================================
//  DELTA E EM COLISOES pp -- COMPARACAO DE METODOS
// ===========================================================================
//
//  OBJETIVO
//  --------
//  Medir  Delta E = (energia inicial dos feixes) - (energia final dos feixes),
//  e relacionar esse Delta E com a multiplicidade e com a energia cinetica.
//  A ideia fisica: a energia que os feixes "perdem" reaparece como particulas
//  novas (massa) e como movimento (energia cinetica).
//
//  O PROBLEMA
//  ----------
//  A energia inicial e' trivial (14000 GeV). A energia FINAL dos feixes nao e',
//  porque em HardQCD os protons NAO sobrevivem intactos: eles se fragmentam.
//  Verificamos isso listando as particulas finais mais energeticas -- os protons
//  de maior energia tinham ~3000 GeV (nao os 7000 iniciais) e rapidez y ~ +-8.7.
//  Ou seja: o proton se quebra, mas os fragmentos continuam na direcao do feixe.
//  Logo, "energia final do feixe" precisa de uma DEFINICAO. Testamos quatro.
//
//  OS QUATRO METODOS
//  -----------------
//  1) BARION LIDER (cinematico) -- FUNCIONA
//     Pega o proton/neutron mais energetico indo para frente e o mais energetico
//     indo para tras. Justificativa: o Pythia constroi o remanescente do feixe
//     conservando o numero barionico do proton, entao esse numero segue adiante
//     carregado por um barion lider. Nao depende de nenhum corte arbitrario.
//
//  2 e 3) CORTE DE RAPIDEZ (cinematico) -- FUNCIONA, mas depende do corte
//     Soma a energia de TODAS as finais com |y| acima de um corte. Justificativa:
//     quem continua na direcao do feixe herda rapidez alta; quem foi criado na
//     colisao sai perpendicular, com y perto de zero.
//     LIMITACAO MEDIDA: o resultado e' muito sensivel ao corte, porque a energia
//     cresce exponencialmente com y (E = mT*cosh(y)). Medimos:
//        |y|>6 -> 1004 GeV | |y|>7 -> 2273 | |y|>8 -> 5378 | |y|>9 -> 10950
//
//  4) STATUS CODE / GENEALOGIA -- NAO FUNCIONA (e este arquivo mostra por que)
//     Ideia: o Pythia etiqueta o remanescente do feixe com status 63. Entao,
//     para cada particula final, subir a arvore de maes procurando um ancestral
//     com essa etiqueta.
//     OBSTACULO 1 (resolvido): a TParticle do ROOT converte o status para 0/1
//        (nao-final / final) e descarta o codigo original. Solucao: acessar o
//        registro nativo com pythia8->Pythia8(). Os paineis 1 e 2 do Canvas 2
//        mostram a mesma informacao pelas duas interfaces.
//     OBSTACULO 2 (fatal): medimos que ~59% das particulas finais descendem
//        do remanescente E do processo duro AO MESMO TEMPO. Motivo fisico: na
//        hadronizacao, a string de cor liga o remanescente aos partons do
//        espalhamento duro, e os hadrons nascem ao longo dela. Por isso a
//        "mae" de um hadron e' o extremo da string, nao um pai fisico -- tanto
//        que as cadeias chegam a cruzar entre os dois feixes.
//     CONCLUSAO: em QCD, atribuir um hadron final a um parton inicial nao e'
//        uma pergunta bem posta. Os metodos cinematicos sao a via operacional.
//
// ===========================================================================

#include "TSystem.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "Pythia8/Pythia.h"   // acesso ao registro nativo (status codes reais)
#include <vector>

// Ajusta a faixa VISIVEL dos eixos em torno de onde os dados estao.
// Nao altera o conteudo do histograma, so o zoom do desenho.
void zoom(TH2F *h)
{
  Double_t mx = h->GetMean(1), sx = h->GetRMS(1);
  Double_t my = h->GetMean(2), sy = h->GetRMS(2);
  h->GetXaxis()->SetRangeUser(TMath::Max(h->GetXaxis()->GetXmin(), mx - 4*sx),
                              TMath::Min(h->GetXaxis()->GetXmax(), mx + 4*sx));
  h->GetYaxis()->SetRangeUser(TMath::Max(h->GetYaxis()->GetXmin(), my - 4*sy),
                              TMath::Min(h->GetYaxis()->GetXmax(), my + 4*sy));
}

// nev        = eventos para os metodos cinematicos (rapidos)
// nev_status = eventos para a analise de genealogia (lenta: percorre a arvore
//              de ancestrais de cada particula final)
void teste_investiga_protons(Int_t nev = 10000, Int_t nev_status = 500)
{
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  Double_t E_inicial = 14000;   // energia total dos dois feixes (GeV)
  Double_t corte_y1  = 8.0;
  Double_t corte_y2  = 9.0;

  // ---------- METODO 1: barion lider ----------
  TH2F *hist_mult = new TH2F("hist_mult","1) Multiplicidade vs #DeltaE (barion lider);#DeltaE (GeV);N particulas",100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec   = new TH2F("hist_ec","1) Energia cinetica vs #DeltaE (barion lider);#DeltaE (GeV);Soma Ecin (GeV)",100, 0, 14000, 100, 13800, 14000);

  // ---------- METODO 2: |y| > 8 ----------
  TH2F *hist_mult_y1 = new TH2F("hist_mult_y1","2) Multiplicidade vs #DeltaE (|y|>8);#DeltaE (GeV);N particulas",100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec_y1   = new TH2F("hist_ec_y1","2) Energia cinetica vs #DeltaE (|y|>8);#DeltaE (GeV);Soma Ecin (GeV)",100, 0, 14000, 100, 13800, 14000);

  // ---------- METODO 3: |y| > 9 ----------
  TH2F *hist_mult_y2 = new TH2F("hist_mult_y2","3) Multiplicidade vs #DeltaE (|y|>9);#DeltaE (GeV);N particulas",100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec_y2   = new TH2F("hist_ec_y2","3) Energia cinetica vs #DeltaE (|y|>9);#DeltaE (GeV);Soma Ecin (GeV)",100, 0, 14000, 100, 13800, 14000);

  // ---------- METODO 4: status 63 + genealogia ----------
  TH2F *hist_mult_st = new TH2F("hist_mult_st","4) Multiplicidade vs #DeltaE (status 63);#DeltaE (GeV);N particulas",100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec_st   = new TH2F("hist_ec_st","4) Energia cinetica vs #DeltaE (status 63);#DeltaE (GeV);Soma Ecin (GeV)",100, 0, 14000, 100, 13800, 14000);

  // ---------- Diagnosticos do metodo 4 ----------
  // Painel 1 x painel 2: a MESMA informacao vista pelas duas interfaces.
  TH1F *hist_status      = new TH1F("hist_status","Status codes REAIS (registro do Pythia);status;particulas", 201, -100.5, 100.5);
  TH1F *hist_status_root = new TH1F("hist_status_root","Status que a TParticle do ROOT entrega;status;particulas", 21, -10.5, 10.5);

  // Classificacao de origem: e' aqui que o metodo 4 se prova inviavel.
  TH1F *hist_origem = new TH1F("hist_origem","De onde descende cada particula final;;fracao", 4, 0, 4);
  hist_origem->GetXaxis()->SetBinLabel(1, "so remanescente");
  hist_origem->GetXaxis()->SetBinLabel(2, "so processo duro");
  hist_origem->GetXaxis()->SetBinLabel(3, "os dois");
  hist_origem->GetXaxis()->SetBinLabel(4, "nenhum");

  TH1F *hist_fracao = new TH1F("hist_fracao","Fracao das finais classificadas como 'do remanescente';fracao;eventos", 102, -0.01, 1.01);

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 43");
  pythia8->Initialize(2212, 2212, 14000);

  // Ponteiro para o objeto Pythia interno. E' o que da acesso aos status codes
  // de verdade -- a TParticle do ROOT nao os repassa.
  Pythia8::Pythia *py = pythia8->Pythia8();

  TClonesArray *particles = new TClonesArray("TParticle", 10000);

  Long64_t total_finais = 0;
  Long64_t cat[4] = {0, 0, 0, 0};   // so remanescente / so duro / os dois / nenhum
  Double_t soma_fracao = 0;

  for (Int_t iev = 0; iev < nev; iev++)
  {
    Double_t soma_ec = 0;
    Int_t multiplicidade = 0;
    Double_t E_lead_fwd = 0, E_lead_bwd = 0;
    Double_t E_feixe_y1 = 0, E_feixe_y2 = 0;

    pythia8->GenerateEvent();

    // =====================================================================
    // LOOP A -- via TParticle (ROOT). Metodos 1, 2 e 3.
    // =====================================================================
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle *)particles->At(ip);

      // Registra o que o ROOT entrega como status (so aparecem 0 e 1).
      if (iev < nev_status) hist_status_root->Fill(part->GetStatusCode());

      // Filtro de estado final: status positivo = sobreviveu ate o fim.
      if (part->GetStatusCode() <= 0) continue;

      Double_t E   = part->Energy();
      Double_t m   = part->GetMass();
      Int_t    pdg = part->GetPdgCode();
      Double_t pz  = part->Pz();
      Double_t y   = part->Y();

      multiplicidade += 1;
      soma_ec += E - m;          // energia cinetica = total - massa de repouso

      // METODO 1: guarda o barion mais energetico de cada sentido.
      // Aceita neutron tambem: em varios eventos o lider de um lado e' neutron,
      // e exigir so proton deixava aquele lado com energia zero (Delta E inflado).
      if (pdg == 2212 || pdg == 2112) {
        if (pz > 0 && E > E_lead_fwd) E_lead_fwd = E;
        if (pz < 0 && E > E_lead_bwd) E_lead_bwd = E;
      }

      // METODOS 2 e 3: soma a energia de tudo com rapidez alta.
      // O teste de Finite descarta os infinitos que aparecem quando E ~ pz
      // (particulas sem massa indo quase exatamente na direcao do feixe).
      if (TMath::Finite(y)) {
        if (TMath::Abs(y) > corte_y1) E_feixe_y1 += E;
        if (TMath::Abs(y) > corte_y2) E_feixe_y2 += E;
      }
    }

    Double_t deltaE    = E_inicial - (E_lead_fwd + E_lead_bwd);
    Double_t deltaE_y1 = E_inicial - E_feixe_y1;
    Double_t deltaE_y2 = E_inicial - E_feixe_y2;

    // So conta o metodo 1 se achou barion lider nos DOIS sentidos, senao o
    // lado sem lider entra com energia zero e o Delta E vai artificialmente
    // para 14000 (o que significaria energia final nula: impossivel).
    if (E_lead_fwd > 0 && E_lead_bwd > 0) {
      hist_mult->Fill(deltaE, multiplicidade);
      hist_ec->Fill(deltaE, soma_ec);
    }

    hist_mult_y1->Fill(deltaE_y1, multiplicidade);
    hist_ec_y1->Fill(deltaE_y1, soma_ec);
    hist_mult_y2->Fill(deltaE_y2, multiplicidade);
    hist_ec_y2->Fill(deltaE_y2, soma_ec);

    // =====================================================================
    // LOOP B -- via registro NATIVO do Pythia. Metodo 4 + diagnosticos.
    // So nos primeiros nev_status eventos, porque percorrer a arvore de
    // ancestrais de cada particula final e' lento.
    // =====================================================================
    if (iev < nev_status)
    {
      Double_t E_remanescente = 0;
      Int_t n_finais = 0, n_do_remanescente = 0;

      for (Int_t i = 0; i < py->event.size(); i++)
      {
        hist_status->Fill(py->event[i].status());   // status REAL

        if (!py->event[i].isFinal()) continue;
        n_finais++;
        total_finais++;

        // ---- sobe a arvore de ancestrais seguindo AS DUAS maes ----
        // Marca se encontrou o remanescente do feixe (63) e/ou a saida do
        // espalhamento duro (23).
        bool tem63 = false, tem23 = false;

        std::vector<int> fila;
        fila.push_back(i);
        Int_t visitados = 0;

        while (!fila.empty() && visitados < 300)
        {
          Int_t idx = fila.back();
          fila.pop_back();
          visitados++;

          if (idx <= 0 || idx >= py->event.size()) continue;

          Int_t st = TMath::Abs(py->event[idx].status());
          if (st == 63) tem63 = true;
          if (st == 23) tem23 = true;

          Int_t m1 = py->event[idx].mother1();
          Int_t m2 = py->event[idx].mother2();

          // Na hadronizacao, m1 e m2 delimitam a STRING inteira (uma faixa de
          // indices), nao duas maes isoladas. E' exatamente por isso que quase
          // todo hadron acaba tendo o remanescente entre os ancestrais.
          if (m2 > m1 && m1 > 0) {
            for (Int_t k = m1; k <= m2 && k - m1 < 50; k++) fila.push_back(k);
          } else {
            if (m1 > 0) fila.push_back(m1);
            if (m2 > 0 && m2 != m1) fila.push_back(m2);
          }
        }

        if      ( tem63 && !tem23) cat[0]++;
        else if (!tem63 &&  tem23) cat[1]++;
        else if ( tem63 &&  tem23) cat[2]++;
        else                       cat[3]++;

        if (tem63) {
          E_remanescente += py->event[i].e();
          n_do_remanescente++;
        }
      }

      Double_t deltaE_st = E_inicial - E_remanescente;
      Double_t fracao = (n_finais > 0) ? (Double_t)n_do_remanescente / n_finais : 0;

      soma_fracao += fracao;
      hist_fracao->Fill(fracao);
      hist_mult_st->Fill(deltaE_st, multiplicidade);
      hist_ec_st->Fill(deltaE_st, soma_ec);

      if (iev < 10)
        printf("Evento %2d:  barion=%8.1f  |y|>8=%8.1f  |y|>9=%8.1f  status63=%8.1f  mult=%3d\n",
               iev, deltaE, deltaE_y1, deltaE_y2, deltaE_st, multiplicidade);
    }
  }

  // =======================================================================
  // RESUMO
  // =======================================================================
  printf("\n=================== DELTA E MEDIO ===================\n");
  printf("  1) barion lider : %8.1f GeV  (RMS %7.1f)  <- funciona\n", hist_mult->GetMean(1),    hist_mult->GetRMS(1));
  printf("  2) |y| > %.1f     : %8.1f GeV  (RMS %7.1f)  <- funciona\n", corte_y1, hist_mult_y1->GetMean(1), hist_mult_y1->GetRMS(1));
  printf("  3) |y| > %.1f     : %8.1f GeV  (RMS %7.1f)  <- funciona\n", corte_y2, hist_mult_y2->GetMean(1), hist_mult_y2->GetRMS(1));
  printf("  4) status 63    : %8.1f GeV  (RMS %7.1f)  <- FALHA\n",    hist_mult_st->GetMean(1), hist_mult_st->GetRMS(1));
  printf("  (RMS grande em relacao a media = definicao instavel)\n");

  printf("\n============ POR QUE O METODO 4 FALHA ==============\n");
  const char *nomes[4] = {"so remanescente ", "so processo duro", "os dois         ", "nenhum          "};
  for (Int_t k = 0; k < 4; k++) {
    Double_t frac = (total_finais > 0) ? (Double_t)cat[k] / total_finais : 0;
    printf("  %s : %.1f%%\n", nomes[k], 100 * frac);
    hist_origem->SetBinContent(k + 1, frac);
  }
  printf("\n  A maioria das particulas finais descende das DUAS origens ao\n");
  printf("  mesmo tempo, porque a string de cor liga o remanescente do feixe\n");
  printf("  aos partons do espalhamento duro. Logo o criterio nao separa\n");
  printf("  'feixe' de 'produzido' -- a pergunta nao e' bem posta em QCD.\n");
  printf("  Fracao media marcada como 'do remanescente': %.3f\n", soma_fracao / nev_status);
  printf("====================================================\n");

  // =======================================================================
  // CANVAS 1 -- os tres metodos cinematicos (funcionam)
  // =======================================================================
  zoom(hist_mult);    zoom(hist_ec);
  zoom(hist_mult_y1); zoom(hist_ec_y1);
  zoom(hist_mult_y2); zoom(hist_ec_y2);

  TCanvas *c1 = new TCanvas("c1", "Metodos cinematicos", 1600, 1400);
  c1->Divide(2, 3);
  c1->cd(1); hist_mult->Draw("COLZ");
  c1->cd(2); hist_ec->Draw("COLZ");
  c1->cd(3); hist_mult_y1->Draw("COLZ");
  c1->cd(4); hist_ec_y1->Draw("COLZ");
  c1->cd(5); hist_mult_y2->Draw("COLZ");
  c1->cd(6); hist_ec_y2->Draw("COLZ");
  c1->SaveAs("deltaE_cinematico.png");

  // =======================================================================
  // CANVAS 2 -- a investigacao via status code (e por que falha)
  // =======================================================================
  zoom(hist_mult_st); zoom(hist_ec_st);

  TCanvas *c2 = new TCanvas("c2", "Investigacao via status code", 1600, 1400);
  c2->Divide(2, 3);

  // Paineis 1 e 2: mesma informacao, duas interfaces.
  c2->cd(1); gPad->SetLogy(); hist_status->Draw();
  c2->cd(2); gPad->SetLogy(); hist_status_root->Draw();

  // Painel 3: o resultado que inviabiliza o metodo.
  c2->cd(3); hist_origem->SetFillColor(38); hist_origem->SetStats(0); hist_origem->Draw("BAR");

  // Painel 4: em muitos eventos, 100% das finais sao marcadas como "do feixe".
  c2->cd(4); gPad->SetLogy(); hist_fracao->Draw();

  // Paineis 5 e 6: o Delta E resultante, degradado.
  c2->cd(5); hist_mult_st->Draw("COLZ");
  c2->cd(6); hist_ec_st->Draw("COLZ");

  c2->SaveAs("deltaE_status.png");
}