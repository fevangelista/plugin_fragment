#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

using namespace boost;

namespace psi{ namespace plugin_fragment {

void localize_on_fragment(Options& options,boost::shared_ptr<Wavefunction> wfn)
{
    int print = options.get_int("PRINT");


    // Compute the overlap matrix
    boost::shared_ptr<IntegralFactory> integral_ = wfn->integral();

    boost::shared_ptr<BasisSet> basisset_ = wfn->basisset();
    boost::shared_ptr<Molecule> mol = basisset_->molecule();
    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S_ao(new Matrix("S_ao",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S_ao);

//    S_ao->print();
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    SharedMatrix AO2SO_ = pet->aotoso();
//    AO2SO_->print();

    // Compute the W matrix for each fragment
    std::vector<int> flist({0});
    std::vector<int> glist;
    boost::shared_ptr<Molecule> frag = mol->extract_subsets(flist,glist);
    int frag_natom = frag->natom();
    fprintf(outfile,"\n  Fragment contains %d atoms",frag->natom());

    // Form a copy of S_ao and zero the rows and columns that are not on this fragment
//    max_a = min_a +
    int nbf = basisset_->nbf();
    int nbfA = 0;
    SharedMatrix S_f(S_ao->clone());
    for (int mu = 0; mu < nbf; mu++) {
        int A = basisset_->function_to_center(mu);
        fprintf(outfile,"\n  Function %d is on atom %d",mu,A);
        if (A < frag_natom){
            nbfA += 1;
        }
    }
    S_f = SharedMatrix(new Matrix("S_f",nbfA,nbfA));
    for (int mu = 0; mu < nbfA; mu++) {
        for (int nu = 0; nu < nbfA; nu++) {
            S_f->set(mu,nu,S_ao->get(mu,nu));
        }
    }
//    S_f->print();
    SharedMatrix S_f_inv(S_f->clone());
    SharedMatrix id(S_f->clone());
    S_f_inv->general_invert();
    SharedMatrix S_f_inv_large(S_ao->clone());
    S_f_inv_large->zero();
    for (int mu = 0; mu < nbfA; mu++) {
        for (int nu = 0; nu < nbfA; nu++) {
            S_f_inv_large->set(mu,nu,S_f_inv->get(mu,nu));
        }
    }
//    S_f_inv_large->print();
    SharedMatrix Ca = wfn->Ca();
    S_f_inv_large->transform(S_ao);
    S_f_inv_large->transform(Ca);
//    S_f_inv_large->print();

    int nalpha = wfn->nalpha();
    SharedMatrix SAo(new Matrix("SAo",nalpha,nalpha));
    for (int i = 0; i < nalpha; ++i){
        for (int j = 0; j < nalpha; ++j){
            SAo->set(i,j,S_f_inv_large->get(i,j));
        }
    }
    SharedMatrix Uo(new Matrix("Uo",nalpha,nalpha));
    SharedVector lo(new Vector("lo",nalpha));
    SAo->diagonalize(Uo,lo,descending);
    // lo->print();

    int navir = nbf - nalpha;
    SharedMatrix SAv(new Matrix("SAv",navir,navir));
    for (int i = 0; i < navir; ++i){
        for (int j = 0; j < navir; ++j){
            SAv->set(i,j,S_f_inv_large->get(i + nalpha,j + nalpha));
        }
    }
    SharedMatrix Uv(new Matrix("Uv",navir,navir));
    SharedVector lv(new Vector("lv",navir));
    SAv->diagonalize(Uv,lv,descending);
    // lv->print();

    SharedMatrix U(new Matrix("U",nbf,nbf));
    for (int i = 0; i < nalpha; ++i){
        for (int j = 0; j < nalpha; ++j){
            U->set(i,j,Uo->get(i,j));
        }
    }
    for (int a = 0; a < navir; ++a){
        for (int b = 0; b < navir; ++b){
            U->set(a + nalpha,b + nalpha,Uv->get(a,b));
        }
    }
//    U->print();
    SharedMatrix Ca_new(Ca->clone());
    Ca_new->gemm(false, false, 1.0,Ca,U, 0.0);
    Ca->copy(Ca_new);

    S_ao->transform(Ca);
//    S_ao->print();
}

}} // End namespaces




