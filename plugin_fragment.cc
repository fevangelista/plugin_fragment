#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace plugin_fragment {

void localize_on_fragment(Options& options,boost::shared_ptr<Wavefunction> wfn);
void localize_on_atoms(Options& options,boost::shared_ptr<Wavefunction> wfn);

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "PLUGIN_FRAGMENT"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType plugin_fragment(Options& options)
{
    // Get the number of fragments from the wfn object
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Molecule> mol = wfn->molecule();
    int nfragments = mol->nfragments();

    // nfragments = 1
    // localize on atoms

    if (options.get_int("PRINT") > 0){
        outfile->Printf("\n\n         ---------------------------------------------------------");
        outfile->Printf("\n                            ORBITAL LOCALIZATION");
        outfile->Printf("\n                        by Francesco A. Evangelista");
        outfile->Printf("\n         ---------------------------------------------------------");

        outfile->Printf("\n\n  ==> Details <==");
        outfile->Printf("\n\n    Number of fragments: %d",nfragments);
    }

    if (nfragments == 1){
        localize_on_atoms(options,wfn);
    }else{
        localize_on_fragment(options,wfn);
    }
    return Success;
}

}} // End namespaces



