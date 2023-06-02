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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "descriptor/bcl_descriptor_atom_is_sp3.h"
#include "descriptor/bcl_descriptor_atom_relative_property_score.h"
#include "descriptor/bcl_descriptor_molecule_similarity.h"
#include "descriptor/bcl_descriptor_pair_convolution_correlation_dnn.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_cheminfo_properties.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_atom_aromaticity_axes.h"
#include "descriptor/bcl_descriptor_atom_effective_polarizability.h"
#include "descriptor/bcl_descriptor_atom_estimated_surface_area.h"
#include "descriptor/bcl_descriptor_atom_formal_charge.h"
#include "descriptor/bcl_descriptor_atom_hbond_info.h"
#include "descriptor/bcl_descriptor_atom_is_sp.h"
#include "descriptor/bcl_descriptor_atom_is_sp2.h"
#include "descriptor/bcl_descriptor_atom_lone_pair_en.h"
#include "descriptor/bcl_descriptor_atom_neighbor_direction.h"
#include "descriptor/bcl_descriptor_atom_number_valences.h"
#include "descriptor/bcl_descriptor_atom_pi_charge.h"
#include "descriptor/bcl_descriptor_atom_polarizability.h"
#include "descriptor/bcl_descriptor_atom_ring_size.h"
#include "descriptor/bcl_descriptor_atom_stereocenters.h"
#include "descriptor/bcl_descriptor_atom_surface_area.h"
#include "descriptor/bcl_descriptor_atom_topological_polar_surface_area.h"
#include "descriptor/bcl_descriptor_atom_type_number.h"
#include "descriptor/bcl_descriptor_atom_type_property_retriever.h"
#include "descriptor/bcl_descriptor_atom_vcharge.h"
#include "descriptor/bcl_descriptor_atom_volume.h"
#include "descriptor/bcl_descriptor_atomic_number.h"
#include "descriptor/bcl_descriptor_binary_operation.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "descriptor/bcl_descriptor_element_type_property_retriever.h"
#include "descriptor/bcl_descriptor_molecule_atom_environment_map.h"
#include "descriptor/bcl_descriptor_molecule_complexity.h"
#include "descriptor/bcl_descriptor_molecule_druglike.h"
#include "descriptor/bcl_descriptor_molecule_entropy_qha.h"
#include "descriptor/bcl_descriptor_molecule_girth.h"
#include "descriptor/bcl_descriptor_molecule_halogenated_aromatic_rings.h"
#include "descriptor/bcl_descriptor_molecule_lipinski_violations.h"
#include "descriptor/bcl_descriptor_molecule_log_p.h"
#include "descriptor/bcl_descriptor_molecule_log_p2008.h"
#include "descriptor/bcl_descriptor_molecule_one_four_clash_score.h"
#include "descriptor/bcl_descriptor_molecule_rings.h"
#include "descriptor/bcl_descriptor_molecule_rotatable_bonds.h"
#include "descriptor/bcl_descriptor_molecule_total_bond_energy.h"
#include "descriptor/bcl_descriptor_molecule_vdw_score.h"
#include "descriptor/bcl_descriptor_molecule_xlog_p.h"
#include "descriptor/bcl_descriptor_molecule_xpka.h"
#include "descriptor/bcl_descriptor_named.h"
#include "descriptor/bcl_descriptor_protein_ligand_correlation_dnn.h"
#include "descriptor/bcl_descriptor_sequence_size.h"
#include "descriptor/bcl_descriptor_sequence_statistics.h"
#include "math/bcl_math_assignment_by_comparison.h"
#include "math/bcl_math_plus_equals.h"
#include "math/bcl_math_running_min_max.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // A single instance of the class, ensures this class is instantiated
    const CheminfoProperties &CheminfoProperties::s_Instance( GetCheminfoProperties());

    //! construct all CheminfoProperties
    CheminfoProperties::CheminfoProperties() :
      calc_SigmaValenceStateIonizationPotential
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_SigmaValenceStateIonizationPotential))
      ),
      calc_SigmaValenceStateElectronAffinity
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_SigmaValenceStateElectronAffinity))
      ),
      calc_SigmaOrbitalElectronegativityMulliken
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_SigmaOrbitalElectronegativityMulliken))
      ),
      calc_SigmaOrbitalElectronegativityPauling
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_SigmaOrbitalElectronegativityPauling))
      ),
      calc_PiValenceStateIonizationPotential
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_PiValenceStateIonizationPotential))
      ),
      calc_PiValenceStateElectronAffinity
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_PiValenceStateElectronAffinity))
      ),
      calc_PiOrbitalElectronegativityMulliken
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_PiOrbitalElectronegativityMulliken))
      ),
      calc_PiOrbitalElectronegativityPauling
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_PiOrbitalElectronegativityPauling))
      ),
      calc_LonePairIonizationPotential
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_LonePairIonizationPotential))
      ),
      calc_LonePairElectronAffinity
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_LonePairElectronAffinity))
      ),
      calc_LonePairElectronegativityMulliken
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_LonePairElectronegativityMulliken))
      ),
      calc_AdditiveAtomicPolarizability
      (
        AddEnum( new AtomTypePropertyRetriever( chemistry::AtomTypeData::e_AdditiveAtomicPolarizability))
      ),
      calc_Polarizability( AddEnum( new AtomPolarizability())),
      calc_FormalCharge( AddEnum( new AtomFormalCharge())),
      calc_AtomMinRingSize( AddEnum( new AtomRingSize( false))),
      calc_AtomMaxRingSize( AddEnum( new AtomRingSize( true))),
      calc_AtomIsSP3( AddEnum( new AtomIsSP3)),
      calc_AtomIsSP2( AddEnum( new AtomIsSP2)),
      calc_AtomIsSP( AddEnum( new AtomIsSP)),
      calc_Mass( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_Mass))),
      calc_GyromagneticRatio( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_GyromagneticRatio))),
      calc_CovalentRadius( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_CovalentRadius))),
      calc_VDWaalsRadius( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_VDWaalsRadius))),
      calc_MeltingPoint( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_MeltingPoint))),
      calc_BoilingPoint( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_BoilingPoint))),
      calc_ElectroNegativity( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_ElectroNegativity))),
      calc_IonizationPotential( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_IonizationPotential))),
      calc_MainGroup( AddEnum( new ElementTypePropertyRetriever( chemistry::ElementTypeData::e_MainGroup))),
      calc_HbondAcceptor( AddEnum( new AtomHBondInfo( AtomHBondInfo::e_Acceptor))),
      calc_HbondDonor( AddEnum( new AtomHBondInfo( AtomHBondInfo::e_Donor))),
      calc_Identity
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            Constants< chemistry::AtomConformationalInterface, float>( 1.0),
            "Atom_Identity",
            "1, for any atom or molecule. The Atom_ prefix is purely for backwards-compatibility. "
            "In new descriptor files, the use of Constant(1) is preferred"
          )
        )
      ),
      calc_NumberValences
      (
        AddEnum
        (
          new AtomNumberValences()
        )
      ),
      calc_AtomicNumber( AddEnum( new AtomicNumber())),
      calc_AtomTypeNumber( AddEnum( new AtomTypeNumber())),
      calc_IsH
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::less>(),
              Constants< chemistry::AtomConformationalInterface, float>( 1.5)
            ),
            "IsH",
            "Returns 1 for hydrogen atoms, 0 for heavy atoms"
          )
        )
      ),
      calc_IsNotH
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::greater>(),
              Constants< chemistry::AtomConformationalInterface, float>( 1.5)
            ),
            "IsNotH",
            "Returns 1 for heavy atoms, 0 for hydrogen atoms"
          )
        )
      ),
      calc_IsC
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 6.0)
            ),
            "IsC",
            "Returns 1 for carbon atoms, 0 for others"
          )
        )
      ),
      calc_IsNotC
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::not_equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 6.0)
            ),
            "IsNotC",
            "Returns 1 for non-carbon atoms, 0 for others"
          )
        )
      ),
      calc_IsB            // 1 for B only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 5.0)
            ),
            "IsB",
            "Returns 1 for boron atoms, 0 for others"
          )
        )
      ),
      calc_IsN            // 1 for N only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 7.0)
            ),
            "IsN",
            "Returns 1 for nitrogen atoms, 0 for others"
          )
        )
      ),
      calc_IsO            // 1 for O only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 8.0)
            ),
            "IsO",
            "Returns 1 for oxygen atoms, 0 for others"
          )
        )
      ),
      calc_IsP            // 1 for O only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 15.0)
            ),
            "IsP",
            "Returns 1 for phosphorus atoms, 0 for others"
          )
        )
      ),
      calc_IsF            // 1 for F only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 9.0)
            ),
            "IsF",
            "Returns 1 for fluorine atoms, 0 for others"
          )
        )
      ),
      calc_IsS            // 1 for S only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 16.0)
            ),
            "IsS",
            "Returns 1 for sulfur atoms, 0 for others"
          )
        )
      ),
      calc_IsCl           // 1 for Cl only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 17.0)
            ),
            "IsCl",
            "Returns 1 for chlorine atoms, 0 for others"
          )
        )
      ),
      calc_IsBr           // 1 for Cl only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 35.0)
            ),
            "IsBr",
            "Returns 1 for bromine atoms, 0 for others"
          )
        )
      ),
      calc_IsI           // 1 for I only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 53.0)
            ),
            "IsI",
            "Returns 1 for iodine atoms, 0 for others"
          )
        )
      ),
      calc_IsSi           // 1 for I only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 14.0)
            ),
            "IsSi",
            "Returns 1 for silicon atoms, 0 for others"
          )
        )
      ),
      calc_IsHalogen      // 1 for Halogens only
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_MainGroup,
              math::AssignmentByComparison< float, std::equal_to>(),
              Constants< chemistry::AtomConformationalInterface, float>( 7.0)
            ),
            "IsHalogen",
            "Returns 1 for atoms in main group 7 (F,Cl,Br,I,At,Uus)"
          )
        )
      ),
      calc_IsPeriod3Plus // 1 for atoms at or below the third row of the periodic table
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_AtomicNumber,
              math::AssignmentByComparison< float, std::greater>(),
              Constants< chemistry::AtomConformationalInterface, float>( 10.5)
            ),
            "IsPeriodThreePlus",
            "Returns 1 for all elements in periods 3-7"
          )
        )
      ),
      calc_VdwVolume
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            AtomVolume
            (
              calc_VDWaalsRadius,
              calc_CovalentRadius
            ),
            "Atom_VDWVolume",
            "approximates the volume of atoms using the vdw radius, considering overlap from neighboring atoms"
          )
        )
      ),
      calc_CovalentVolume
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            AtomVolume
            (
              calc_CovalentRadius,
              CheminfoProperty()
            ),
            "Atom_CovalentVolume",
            "approximates the volume of atoms using the covalent radius, considering overlap from neighboring atoms"
          )
        )
      ),
      calc_VdwSurfaceArea
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            AtomSurfaceArea
            (
              calc_VDWaalsRadius,
              calc_CovalentRadius
            ),
            "Atom_VDWSurfaceArea",
            "approximates the surface area of atoms using the vdw radius, considering overlap from neighboring atoms"
          )
        )
      ),
      calc_CovalentSurfaceArea
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            AtomSurfaceArea
            (
              calc_CovalentRadius,
              CheminfoProperty()
            ),
            "Atom_CovalentSurfaceArea",
            "approximates the surface area of atoms using the covalent radius, considering overlap from neighboring atoms"
          )
        )
      ),

      calc_NeighborDirection
      (
        AddEnum
        (
          new AtomNeighborDirection()
        )
      ),
      calc_AromaticityAxes
            (
              AddEnum
              (
                new AtomAromaticityAxes()
              )
            ),
      calc_EstVdwSurfaceArea( AddEnum( new AtomEstimatedSurfaceArea( false, false))),
      calc_EstVdwSurfaceAreaCSD( AddEnum( new AtomEstimatedSurfaceArea( false, true))),
      calc_EstCovalentSurfaceArea( AddEnum( new AtomEstimatedSurfaceArea( true))),
      calc_LonePairEN( AddEnum( new AtomLonePairEN())),
      calc_PiCharge( AddEnum( new AtomPiCharge( true))),
      calc_PiEN( AddEnum( new AtomPiCharge( false))),
      calc_SigmaCharge( AddEnum( new AtomSigmaCharge( true))),
      calc_SigmaEN( AddEnum( new AtomSigmaCharge( false))),
      calc_VCharge( AddEnum( new AtomVcharge( 0))),
      calc_VChargeV2( AddEnum( new AtomVcharge( 1))),
      calc_TotalCharge
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_SigmaCharge,
              math::PlusEquals< float>(),
              calc_PiCharge
            ),
            "Atom_TotalCharge",
            "Returns the total charge on an atom"
          )
        )
      ),
      calc_EffectivePolarizability( AddEnum( new AtomEffectivePolarizability())),
      calc_TopologicalPolarSurfaceArea( AddEnum( new AtomTopologicalPolarSurfaceArea())),
      calc_Stereocenters( AddEnum( new AtomStereocenters())),
      calc_Girth( AddEnum( new MoleculeGirth())),
      calc_MolHbondAcceptor
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_HbondAcceptor
            ),
            "HbondAcceptor",
            "# hydrogen bond acceptors"
          )
        )
      ),
      calc_MolHbondDonor
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_HbondDonor
            ),
            "HbondDonor",
            "# hydrogen bond donors"
          )
        )
      ),
      calc_LogP( AddEnum( new MoleculeLogP())),
      calc_LogP2008( AddEnum( new MoleculeLogP2008())),
      calc_XLogP( AddEnum( new MoleculeXLogP)),
      calc_XpKaAcid( AddEnum( new MoleculeXpKa( true))),
      calc_XpKaBase( AddEnum( new MoleculeXpKa( false))),
      calc_NAtoms
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceSize< chemistry::AtomConformationalInterface, float>(),
            "NAtoms",
            "Number of atoms"
          )
        )
      ),
      calc_NHeavyAtoms
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_IsNotH
            ),
            "NHeavyAtoms",
            "Number of non-hydrogen atoms in the molecule"
          )
        )
      ),
      calc_NRotBond( AddEnum( new MoleculeRotatableBonds())),
      calc_NRotBondSym( AddEnum( new MoleculeRotatableBonds( true))),
      calc_NStereo
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              BinaryOperation< chemistry::AtomConformationalInterface>
              (
                AtomStereocenters(),
                math::AssignmentByComparison< float, std::not_equal_to>(),
                Constants< chemistry::AtomConformationalInterface, float>( 0.0)
              )
            ),
            "NStereo",
            "Number of stereocenters"
          )
        )
      ),
      calc_MolPolarizability
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_Polarizability
            ),
            "Polarizability",
            "total polarizability"
          )
        )
      ),
      calc_MolTotalCharge
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_TotalCharge
            ),
            "TotalCharge",
            "Sum of sigma and pi charges"
          )
        )
      ),
      calc_MolTotalFormalCharge
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_FormalCharge
            ),
            "TotalFormalCharge",
            "Sum of formal charges on the molecule"
          )
        )
      ),
      calc_MolWeight
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_Mass
            ),
            "Weight",
            "Molecular weight (amu)"
          )
        )
      ),
      calc_MolComplexity( AddEnum( new MoleculeComplexity)),
      calc_MolLipinskiViolations( AddEnum( new MoleculeLipinskiViolations)),
      calc_MolLipinskiDruglike
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            BinaryOperation< chemistry::AtomConformationalInterface>
            (
              calc_MolLipinskiViolations,
              math::AssignmentByComparison< float, std::less>(),
              Constants< chemistry::AtomConformationalInterface, float>( 2)
            ),
            "LipinskiDruglike",
            "Returns 1 if the number of Lipinski violations is less than 2"
          )
        )
      ),
      calc_MolLipinskiViolationsVeber( AddEnum( new MoleculeLipinskiViolations( MoleculeLipinskiViolations::e_Veber))),
      calc_MolAromaticRingHalogensTotal( AddEnum( new MoleculeHalogenatedAromaticRings( MoleculeHalogenatedAromaticRings::e_Total))),
      calc_MolAromaticRingHalogensMax( AddEnum( new MoleculeHalogenatedAromaticRings( MoleculeHalogenatedAromaticRings::e_Max))),
      calc_MolEntropyQHA( AddEnum( new MoleculeEntropyQHA)),
      calc_MolTotalBondEnergy( AddEnum( new MoleculeTotalBondEnergy)),
      calc_IsMolDruglike( AddEnum( new MoleculeDruglike( false))),
      calc_IsMolDruglikeAndHitlike( AddEnum( new MoleculeDruglike( true))),
      calc_AffinityNet( AddEnum( new ProteinLigandCorrelationDNN( false, false))),
      calc_AffinityNetAD( AddEnum( new ProteinLigandCorrelationDNN( true, false))),
      calc_DockANNScore( AddEnum( new ProteinLigandCorrelationDNN( false, true))),
      calc_PCCBindingAffinity( AddEnum( new PairConvolutionCorrelationDNN( false))),
      calc_PCCBindingAffinityWeighted( AddEnum( new PairConvolutionCorrelationDNN( true))),
      calc_MoleculeVdwScore( AddEnum( new MoleculeVdwScore)),
      calc_MoleculeOneFourClashScore( AddEnum( new MoleculeOneFourClashScore)),
      calc_MoleculeAtomEnvironmentMap( AddEnum( new MoleculeAtomEnvironmentMap)),
      calc_AtomRelativePropertyScore( AddEnum( new AtomRelativePropertyScore)),
      calc_MolCovalentSurfaceArea
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_CovalentSurfaceArea
            ),
            "CovalentSurfaceArea",
            "covalent surface area"
          )
        )
      ),
      calc_MolCovalentVolume
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_CovalentVolume
            ),
            "CovalentVolume",
            "covalent volume"
          )
        )
      ),
      calc_MolTopologicalPolarSurfaceArea
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_TopologicalPolarSurfaceArea
            ),
            "TopologicalPolarSurfaceArea",
            "topological polar surface area"
          )
        )
      ),
      calc_MolVdwSurfaceArea
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_VdwSurfaceArea
            ),
            "VdwSurfaceArea",
            "van der waals surface area"
          )
        )
      ),
      calc_MolVdwVolume
      (

        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_VdwVolume
            ),
            "VdwVolume",
            "van der waals volume"
          )
        )
      ),
      calc_MolEstCovalentSurfaceArea
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_EstCovalentSurfaceArea
            ),
            "EstCovSurfaceArea",
            "Estimated covalent surface area (conformation independent)"
          )
        )
      ),
      calc_MolEstVdwSurfaceArea
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_EstVdwSurfaceArea
            ),
            "EstVdwSurfaceArea",
            "Estimated van-der-waals surface area (conformation independent), using element-based VDW radii"
          )
        )
      ),
      calc_MolEstVdwSurfaceAreaCSD
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningSum< linal::Vector< float> >,
              &math::RunningSum< linal::Vector< float> >::GetSum
            >
            (
              calc_EstVdwSurfaceAreaCSD
            ),
            "EstVdwSurfaceAreaCSD",
            "Estimated van-der-waals surface area (conformation independent), "
            "using more accurate, CSD-derived atom-type VdW radii, which tend to give larger SAs to H than "
            "EstVdwSurfaceArea"
          )
        )
      ),
      calc_NRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Any))),
      calc_NAromaticRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Aromatic))),
      calc_NConjugatedRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Conjugated))),
      calc_NNonConjugatedRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Nonconjugated))),
      calc_NMacrocyclicRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Any, true))),
      calc_NAromaticMacrocyclicRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Aromatic, true))),
      calc_NConjugatedMacrocyclicRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Conjugated, true))),
      calc_NNonConjugatedMacrocyclicRings( AddEnum( new MoleculeRings( chemistry::ConstitutionalBondTypeData::e_Nonconjugated, true))),
      calc_MinRingSize
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningMinMax< linal::Vector< float> >,
              &math::RunningMinMax< linal::Vector< float> >::GetMin
            >
            (
              calc_AtomMinRingSize
            ),
            "MinRingSize",
            "Smallest ring inside the molecule; returns "
            + util::Format()( size_t( AtomRingSize::s_ChainAtomsMinRingSize)) + " for ringless molecules"
          )
        )
      ),
      calc_MaxRingSize
      (
        AddEnum
        (
          new Named< chemistry::AtomConformationalInterface, float>
          (
            SequenceStatistics
            <
              chemistry::AtomConformationalInterface,
              math::RunningMinMax< linal::Vector< float> >,
              &math::RunningMinMax< linal::Vector< float> >::GetMax
            >
            (
              calc_AtomMaxRingSize
            ),
            "MaxRingSize",
            "Largest unbridged ring inside the molecule; returns "
            + util::Format()( size_t( AtomRingSize::s_ChainAtomsMaxRingSize)) + " for ringless molecules"
          )
        )
      )
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief adds the given enum to the Enumerated instances
    //! @param PROPERTY_PTR pointer to a new property
    //! @return the resulting atom property
    CheminfoProperty CheminfoProperties::AddEnum( Base< chemistry::AtomConformationalInterface, float> *PROPERTY_PTR)
    {
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance( PROPERTY_PTR);
      return CheminfoProperty( PROPERTY_PTR->GetLabel());
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the only instance of atom properties
    //! Overrides GetEnums in enumerate class, which does not return a const-reference
    const CheminfoProperties &CheminfoProperties::GetEnums()
    {
      static CheminfoProperties s_CheminfoProperties;
      return s_CheminfoProperties;
    }

    const CheminfoProperties &GetCheminfoProperties()
    {
      return CheminfoProperties::GetEnums();
    }

  } // namespace descriptor
} // namespace bcl
