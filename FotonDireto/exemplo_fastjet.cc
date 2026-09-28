// Exemplo minimo: Pythia8 gera eventos, FastJet clusteriza em jatos anti-kT.
// Compilar:
//   g++ exemplo_fastjet.cc -o exemplo_fastjet -std=c++17 \
//       -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8 \
//       $(fastjet-config --cxxflags --libs)

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"   // conversao Pythia8::Particle -> fastjet::PseudoJet
#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace Pythia8;

int main() {

  // ---- Pythia ----
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 100.");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 12345");
  if (!pythia.init()) return 1;

  // ---- FastJet: anti-kT, R = 0.4 ----
  double R = 0.4, pTmin = 20.0, etaMax = 4.0;
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);

  int nEvent = 5;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // particulas finais visiveis viram PseudoJets
    std::vector<fastjet::PseudoJet> entradas;
    for (int i = 0; i < pythia.event.size(); ++i) {
      const Particle& p = pythia.event[i];
      if (!p.isFinal()) continue;
      if (!p.isVisible()) continue;              // remove neutrinos
      if (std::abs(p.eta()) > etaMax) continue;
      fastjet::PseudoJet pj(p.px(), p.py(), p.pz(), p.e());
      pj.set_user_index(i);                      // guarda o indice no event record
      entradas.push_back(pj);
    }

    fastjet::ClusterSequence cs(entradas, jetDef);
    std::vector<fastjet::PseudoJet> jatos =
        sorted_by_pt(cs.inclusive_jets(pTmin));

    std::cout << "\n=== Evento " << iEvent << " : " << entradas.size()
              << " particulas -> " << jatos.size() << " jatos (pT > "
              << pTmin << " GeV) ===\n";
    for (size_t j = 0; j < jatos.size(); ++j) {
      std::cout << "  jato " << j
                << "  pT = "  << jatos[j].pt()
                << "  eta = " << jatos[j].eta()
                << "  phi = " << jatos[j].phi()
                << "  m = "   << jatos[j].m()
                << "  n_const = " << jatos[j].constituents().size() << "\n";
    }
  }

  pythia.stat();
  return 0;
}
