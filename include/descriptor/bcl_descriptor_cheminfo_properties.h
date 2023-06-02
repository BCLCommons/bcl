// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_DESCRIPTOR_CHEMINFO_PROPERTIES_H
#define BCL_DESCRIPTOR_CHEMINFO_PROPERTIES_H

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CheminfoProperties
    //! @brief This class enumerates CheminfoProperty
    //!
    //! @see @link example_descriptor_cheminfo_properties.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Feb 14, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CheminfoProperties
    {

    public:

      //! @note the example will output all the names and aliases of all atom properties that currently exist

    //////////////////////////
    // atom type properties //
    //////////////////////////

      CheminfoProperty calc_SigmaValenceStateIonizationPotential; //!< Requires defined atom types
      CheminfoProperty calc_SigmaValenceStateElectronAffinity;    //!< Requires defined atom types
      CheminfoProperty calc_SigmaOrbitalElectronegativityMulliken;//!< Requires defined atom types
      CheminfoProperty calc_SigmaOrbitalElectronegativityPauling; //!< Requires defined atom types
      CheminfoProperty calc_PiValenceStateIonizationPotential;    //!< Requires defined atom types
      CheminfoProperty calc_PiValenceStateElectronAffinity;       //!< Requires defined atom types
      CheminfoProperty calc_PiOrbitalElectronegativityMulliken;   //!< Requires defined atom types
      CheminfoProperty calc_PiOrbitalElectronegativityPauling;    //!< Requires defined atom types
      CheminfoProperty calc_LonePairIonizationPotential;          //!< Requires defined atom types
      CheminfoProperty calc_LonePairElectronAffinity;             //!< Requires defined atom types
      CheminfoProperty calc_LonePairElectronegativityMulliken;    //!< Requires defined atom types
      CheminfoProperty calc_AdditiveAtomicPolarizability;
      CheminfoProperty calc_Polarizability;                       //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_FormalCharge;
      CheminfoProperty calc_AtomMinRingSize;
      CheminfoProperty calc_AtomMaxRingSize;
      CheminfoProperty calc_AtomIsSP3;                            //!< Requires defined atom types
      CheminfoProperty calc_AtomIsSP2;                            //!< Requires defined atom types
      CheminfoProperty calc_AtomIsSP;                            //!< Requires defined atom types

    /////////////////////////////
    // element type properties //
    /////////////////////////////

      CheminfoProperty calc_Mass;
      CheminfoProperty calc_GyromagneticRatio;
      CheminfoProperty calc_CovalentRadius;
      CheminfoProperty calc_VDWaalsRadius;
      CheminfoProperty calc_MeltingPoint;
      CheminfoProperty calc_BoilingPoint;
      CheminfoProperty calc_ElectroNegativity;
      CheminfoProperty calc_IonizationPotential;
      CheminfoProperty calc_MainGroup;

    //////////////////////
    // hydrogen bonding //
    //////////////////////

      CheminfoProperty calc_HbondAcceptor;  //!< 1.0 for hydrogen bond acceptor atoms
      CheminfoProperty calc_HbondDonor;     //!< 1.0 for hydrogen bond donor atoms.  Requires explicit hydrogens

    ////////////////
    // identities //
    ////////////////

      CheminfoProperty calc_Identity;       //!< 1.0 for all atoms
      CheminfoProperty calc_NumberValences; //!< # of valences
      CheminfoProperty calc_AtomicNumber;
      CheminfoProperty calc_AtomTypeNumber;             //!< Index of each AtomType
      CheminfoProperty calc_IsH;                        //!< 1 for H only
      CheminfoProperty calc_IsNotH;                     //!< 1 for non-H only
      CheminfoProperty calc_IsC;            //!< 1 for C only
      CheminfoProperty calc_IsNotC;         //!< 1 for non-C only
      CheminfoProperty calc_IsB;            //!< 1 for B only
      CheminfoProperty calc_IsN;            //!< 1 for N only
      CheminfoProperty calc_IsO;            //!< 1 for O only
      CheminfoProperty calc_IsP;            //!< 1 for P only
      CheminfoProperty calc_IsF;            //!< 1 for F only
      CheminfoProperty calc_IsS;            //!< 1 for S only
      CheminfoProperty calc_IsCl;           //!< 1 for Cl only
      CheminfoProperty calc_IsBr;           //!< 1 for Br only
      CheminfoProperty calc_IsI;            //!< 1 for I only
      CheminfoProperty calc_IsSi;            //!< 1 for C only
      CheminfoProperty calc_IsHalogen;      //!< 1 for Halogens only
      CheminfoProperty calc_IsPeriod3Plus;  //!< 1 for atoms at or below the third row of the periodic table

    ///////////////
    // geometric //
    ///////////////

      CheminfoProperty calc_VdwVolume;
      CheminfoProperty calc_CovalentVolume;
      CheminfoProperty calc_VdwSurfaceArea;
      CheminfoProperty calc_CovalentSurfaceArea;
      CheminfoProperty calc_EstVdwSurfaceArea;
      CheminfoProperty calc_EstVdwSurfaceAreaCSD;
      CheminfoProperty calc_EstCovalentSurfaceArea;
      CheminfoProperty calc_AromaticityAxes;
      CheminfoProperty calc_NeighborDirection;

    //////////
    // misc //
    //////////

      CheminfoProperty calc_LonePairEN;                  //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_PiCharge;                    //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_PiEN;                        //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_SigmaCharge;                 //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_SigmaEN;                     //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_VCharge;                     //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_VChargeV2;                   //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_TotalCharge;                 //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_EffectivePolarizability;     //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_TopologicalPolarSurfaceArea; //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_Stereocenters;               //!< Requires explicit hydrogens

    ////////////////////////////////
    // small molecule properties //
    //////////////////////////////

    /////////////
    // General //
    /////////////

      CheminfoProperty calc_Girth;                         //!< maximum distance between two atoms
      CheminfoProperty calc_MolHbondAcceptor;              //!< number of hydrogen bond acceptors
      CheminfoProperty calc_MolHbondDonor;                 //!< number of hydrogen bond donors
      CheminfoProperty calc_LogP;                          //!< fast log p estimate
      CheminfoProperty calc_LogP2008;                      //!< fast log p estimate, revision over LogP
      CheminfoProperty calc_XpKaAcid;                      //!< multi-tasking DNN prediction for pKa trained on acids
      CheminfoProperty calc_XpKaBase;                      //!< multitasking DNN prediction for pKa trained on bases
      CheminfoProperty calc_XLogP;                         //!< multi-tasking DNN prediction for LogP
      CheminfoProperty calc_NAtoms;                        //!< number of atoms
      CheminfoProperty calc_NHeavyAtoms;                   //!< number of heavy atoms
      CheminfoProperty calc_NRotBond;                      //!< number of rotatable bonds
      CheminfoProperty calc_NRotBondSym;                   //!< number of rotatable bonds, considering symmetry
      CheminfoProperty calc_NStereo;                       //!< number of stereocenters
      CheminfoProperty calc_MolPolarizability;             //!< total polarizability
      CheminfoProperty calc_MolTotalCharge;                //!< sigma + pi charge
      CheminfoProperty calc_MolTotalFormalCharge;          //!< sum of atomic GetFormalCharges
      CheminfoProperty calc_MolWeight;                     //!< molecular weight
      CheminfoProperty calc_MolComplexity;                 //!< Requires defined atom types and explicit hydrogens
      CheminfoProperty calc_MolLipinskiViolations;         //!< Requires explicit hydrogens
      CheminfoProperty calc_MolLipinskiDruglike;           //!< Requires explicit hydrogens
      CheminfoProperty calc_MolLipinskiViolationsVeber;    //!< Requires explicit hydrogens
      CheminfoProperty calc_MolAromaticRingHalogensTotal;  //!< Count the total number of aromatic ring halogen atoms
      CheminfoProperty calc_MolAromaticRingHalogensMax;    //!< Count the max number of aromatic ring halogen atoms on a fragment
      CheminfoProperty calc_MolEntropyQHA;                 //!< Conformational entropy of molecule
      CheminfoProperty calc_MolTotalBondEnergy;            //!< Total bond energy of molecule from component statistical bond energies
      CheminfoProperty calc_IsMolDruglike;                 //!< Whether or not the molecule is druglike
      CheminfoProperty calc_IsMolDruglikeAndHitlike;       //!< Whether or not the molecule is druglike and hitlike
      CheminfoProperty calc_AffinityNet;                   //!< Pose-dependent protein-ligand binding affinity in units of pKd
      CheminfoProperty calc_AffinityNetAD;                 //!< Pose-dependent protein-ligand binding affinity in units of pKd, weighted by AD
      CheminfoProperty calc_DockANNScore;                  //!< Pose-dependent protein-ligand docking score in units of pKd
      CheminfoProperty calc_PCCBindingAffinity;            //!< Pose-independent protein-ligand binding affinity in units of pKd
      CheminfoProperty calc_PCCBindingAffinityWeighted;    //!< Pose-independent protein-ligand binding affinity in units of pKd, weighted by AD
      CheminfoProperty calc_MoleculeVdwScore;              //!< Molecule VDW score from chemistry::AtomVDWScore
      CheminfoProperty calc_MoleculeOneFourClashScore;     //!< Molecule 1-4 score from chemistry::AtomOneFourInteractionScore
      CheminfoProperty calc_MoleculeAtomEnvironmentMap;    //!< Whether or not molecule fragment atom environments are within fragment database
      CheminfoProperty calc_AtomRelativePropertyScore;     //!< Compute the per-atom contribution to a property relative to a congeneric ligand

    ///////////////
    // Geometric //
    ///////////////

      CheminfoProperty calc_MolCovalentSurfaceArea;
      CheminfoProperty calc_MolCovalentVolume;
      CheminfoProperty calc_MolTopologicalPolarSurfaceArea;
      CheminfoProperty calc_MolVdwSurfaceArea;
      CheminfoProperty calc_MolVdwVolume;
      CheminfoProperty calc_MolEstCovalentSurfaceArea;
      CheminfoProperty calc_MolEstVdwSurfaceArea;
      CheminfoProperty calc_MolEstVdwSurfaceAreaCSD;

    ///////////
    // Rings //
    ///////////

      CheminfoProperty calc_NRings;                         //!< Total # of rings (geometric face algorithm)
      CheminfoProperty calc_NAromaticRings;                 //!< Total aromatic rings
      CheminfoProperty calc_NConjugatedRings;               //!< Total conjugated rings
      CheminfoProperty calc_NNonConjugatedRings;            //!< Total non-conjugated rings
      CheminfoProperty calc_NMacrocyclicRings;              //!< Total # of macrocyclic rings (geometric face algorithm)
      CheminfoProperty calc_NAromaticMacrocyclicRings;      //!< Total macrocyclic aromatic rings
      CheminfoProperty calc_NConjugatedMacrocyclicRings;    //!< Total macrocyclic conjugated rings
      CheminfoProperty calc_NNonConjugatedMacrocyclicRings; //!< Total macrocyclic non-conjugated rings
      CheminfoProperty calc_MinRingSize;                    //!< Minimum ring size
      CheminfoProperty calc_MaxRingSize;                    //!< Maximum (non-bridged) ring size

    private:

      //! Single instance of the class, ensures the class is instantiated
      static const CheminfoProperties &s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct all atom property classes
      CheminfoProperties();

    ////////////////
    // operations //
    ////////////////

      //! @brief adds the given enum to the Enumerated instances
      //! @param PROPERTY_PTR pointer to a new property
      //! @return the resulting atom property
      CheminfoProperty AddEnum( Base< chemistry::AtomConformationalInterface, float> *PROPERTY_PTR);

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the only instance of atom properties
      //! Overrides GetEnums in enumerate class, which returns a non-const reference
      //! This is to prevent changing CheminfoProperty enums
      static const CheminfoProperties &GetEnums();

    }; // class CheminfoProperties

    BCL_API
    const CheminfoProperties &GetCheminfoProperties();

  } // namespace descriptor
} // namespace bcl

#endif //BCL_DESCRIPTOR_CHEMINFO_PROPERTIES_H
