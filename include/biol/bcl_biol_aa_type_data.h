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

#ifndef BCL_BIOL_AA_TYPE_DATA_H_
#define BCL_BIOL_AA_TYPE_DATA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_atom_types.h"
#include "bcl_biol_chi_angle.h"
#include "chemistry/bcl_chemistry_atom_types.h"
#include "restraint/bcl_restraint_sas_data_parameters.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AATypeData
    //! @brief This is a low level helper class to store amino acid properties
    //! @details This class acts as the storage class for the enumerator AATypes.
    //! Several additional amino acid scales can be found in ProtScale: http://www.expasy.org/tools/protscale.html
    //!
    //! @see @link example_biol_aa_type_data.cpp @endlink
    //! @author meilerj, woetzen, karakam
    //! @date 10.10.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AATypeData :
      public util::ObjectInterface
    {

    /////////////
    // friends //
    /////////////

      friend class AATypes;

    public:

      //! enum  properties for amino acids
      enum PropertyType
      {
        e_NaturalPrevalence,                      //!< based on natural occurrence
        e_StericalParameter,                      //!< sterical parameter (find correct citation in Meiler et. al. JMM)
        e_Polarizability,                         //!< polarizability     (find correct citation in Meiler et. al. JMM)
        e_Volume,                                 //!< volume             (find correct citation in Meiler et. al. JMM)
        e_Hydrophobicity,                         //!< hydrophobicity     (find correct citation in Meiler et. al. JMM)
        e_IsoelectricPoint,                       //!< isoelectric point  (find correct citation in Meiler et. al. JMM)
        e_Mass,                                   //!< mass of amino acid
        e_Charge,                                 //!< indicate if the amino acid side chain is positively charged = 1.0 (e.g. K) or negatively charged = -1.0 (e.g. E)
        // the following documentation and the pK values were retrieved from http://isoelectric.ovh.org/files/isoelectric-point-theory.html
        e_pK_EMBOSS,                              //!< pK values used in EMBOSS http://emboss.sourceforge.net/apps/cvs/iep.html
        e_pK_DTASelect,                           //!< pK values used in DTASelect http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf
        e_pk_Solomon,                             //!< pK values from Solomon, T.W.G. (1998) Fundamentals of Organic Chemistry, 5th edition. Published by Wiley.
        e_pK_Sillero,                             //!< pK values from Sillero A. and Ribeiro,J.M. (1989) Isoelectric points of proteins: theoretical determination. Anal. Biochem., 179, 319–325. and Sillero A, Maldonado A. (2006) Isoelectric point determination of proteins and other macromolecules: oscillating method. Comput Biol Med. 36(2), 157-66. Epub 2005 Jan 1
        e_pK_Rodwell,                             //!< pK values from Rodwell JD. Heterogeneity of component bands in isoelectric focusing patterns. Anal Biochem. 1982 Jan 15;119(2):440-9. or Murray, R.K., Granner, D.K., Rodwell, V.W. (2006) Harper's illustrated Biochemistry. 27th edition. Published by The McGraw-Hill Companies.
        e_pK_Patrickios,                          //!< pK values from Patrickios CS, Yamasaki EN. Polypeptide amino acid composition and isoelectric point. II. Comparison between experiment and theory. Anal Biochem. 1995 Oct 10;231(1):82-91.
        e_pK_Wikipedia,                           //!< pK values from http://en.wikipedia.org/wiki/Isoelectric_point
        e_pK_Lehninger,                           //!< pK values from Nelson, David L.; Cox, Michael M. (2000). Lehninger Principles of Biochemistry (3rd ed.). Worth Publishers. ISBN 1-57259-153-6.
        e_pK_Grimsely,                            //!< pK values from Grimsley GR, Scholtz JM, Pace CN. A summary of the measured pK values of the ionizable groups in folded proteins. Protein Sci. 2009 Jan;18(1):247-51.
        e_pK_Bjellqvist,                          //!< pK values from Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E. Reference points for comparisons of two-dimensional maps of proteins from different human cell types defined in a pH scale where isoelectric points correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.
        e_pK_ProMoST,                             //!< pk values used in ProMoST Halligan BD. ProMoST: a tool for calculating the pI and molecular mass of phosphorylated and modified proteins on two-dimensional gels. Methods Mol Biol. 2009;527:283-98
        e_pK_Bjellqvist_NTerm,
        e_pK_Bjellqvist_CTerm,
        e_pK_ProMoST_NTerm,
        e_pK_ProMoST_CTerm,
        e_pK_Carey_NTerm,                         //!< N-terminal pKas based on Carey, Organic Chemistry (textbook), http://www.mhhe.com/physsci/chemistry/carey5e/Ch27/ch27-1-4-2.html
        e_pK_Carey_CTerm,                         //!< C-terminal pKas based on Carey, Organic Chemistry (textbook), http://www.mhhe.com/physsci/chemistry/carey5e/Ch27/ch27-1-4-2.html
        e_HelixProbability,                       //!< helix probability
        e_StrandProbability,                      //!< strand probability
        e_FreeEnergyHelix,                        //!< free energy for helix derived from helix probability (from membrane protein database statistics) => see Julia
        e_FreeEnergyStrand,                       //!< free energy for strand derived from strand probability (from membrane protein database statistics) => see Julia
        e_FreeEnergyCoil,                         //!< free energy for coil derived from coil probability (from membrane protein database statistics) => see Julia
        e_TransferFreeEnergyWhimleyWhite,         //!< hydrophobicity scale from White & Wimley,           it represents the transfer free energy for the residues for the transfer from water to octanol.
        e_TransferFreeEnergyEngelmanSeitzGoldman, //!< hydrophobicity scale from Engelman, Steitz, Goldman, it represents the transfer free energy for the residues in helices for the transfer from water to oil.
        e_TransferFreeEnergyKyteDoolittle,        //!< hydrophobicity scale from Kyte and Doolittle        (http://www.expasy.org/tools/protscale.html)
        e_TransferFreeEnergyEisenberg,            //!< hydrophobicity scale from Eisenberg, Schwarz, Wall  (http://www.expasy.org/tools/protscale.html)
        e_TransferFreeEnergyHoppWoods,            //!< hydrophobicity scale from Hopp & Woods, 1981 (all citations can be found in Koehler et. al 2008)
        e_TransferFreeEnergyGuy,                  //!< hydrophobicity scale from Guy, 1985
        e_TransferFreeEnergyJanin,                //!< hydrophobicity scale from Janin, 1979
        e_TransferFreeEnergyPuntaMaritan1D,       //!< hydrophobicity scale from Punta & Maritan 2003
        e_TransferFreeEnergyPuntaMaritan3D,       //!< hydrophobicity scale from Punta & Maritan 2003
        e_FreeEnergyCore,                         //!< our hydrophobicity scale => see Julia's manuscript
        e_FreeEnergyTransition,                   //!< our hydrophobicity scale => see Julia's manuscript
        e_FreeEnergySolution,                     //!< our hydrophobicity scale => see Julia's manuscript
        e_FreeEnergyCoreHelix,                    //!< statistically derived probability for an AA to be in the core region of the membrane and being in helix
        e_FreeEnergyTransitionHelix,              //!< statistically derived probability for an AA to be in the transition region of the membrane and being in helix
        e_FreeEnergySolutionHelix,                //!< statistically derived probability for an AA to be in the solution region of the membrane and being in helix
        e_FreeEnergyCoreStrand,                   //!< statistically derived probability for an AA to be in the core region of the membrane and being in strand
        e_FreeEnergyTransitionStrand,             //!< statistically derived probability for an AA to be in the transition region of the membrane and being in strand
        e_FreeEnergySolutionStrand,               //!< statistically derived probability for an AA to be in the solution region of the membrane and being in strand
        e_FreeEnergyCoreCoil,                     //!< statistically derived probability for an AA to be in the core region of the membrane and being in coil
        e_FreeEnergyTransitionCoil,               //!< statistically derived probability for an AA to be in the transition region of the membrane and being in coil
        e_FreeEnergySolutionCoil,                 //!< statistically derived probability for an AA to be in the solution region of the membrane and being in coil
        e_FreeEnergyCoreFacingPore,               //!< statistically derived probability for an AA to be in the core region of the protein and facing a pore
        e_FreeEnergyCoreFacingMembrane,           //!< statistically derived probability for an AA to be in the core region of the protein and facing a membrane
        e_MembraneStrandOrientationHydrophobicity,//!< statistically derived probability for an AA to be in the core region of the protein and facing a membrane
        e_FreeEnergyExtracellularBlastBB,         //!< statistically derived energy for an AA to be extracellular based on its blast profile (derived from beta-barrels only)
        e_FreeEnergyExtracellularTypeBB,          //!< statistically derived energy for an AA to be extracellular based on its type (derived from beta-barrels only)
        e_SASA,                                   //!< (average of five AAs) solvent accessible surface area calculated by VMD with radius = 1.4 A
        e_SideChainGirth,                         //!< Longest distance between any two side chain atoms, averaged over full atom models of all proteins in jufo training set
        e_SideChainPolarizability,                //!< Polarizability of the side chain using the method described in J. Am. Chem. Soc. 1990, 112, 8533-8542
        e_SideChainTopologicalPolarSurfaceArea,   //!< Topological polar surface area of the side chain using the method from J. Med. Chem. 2000, 43, 3715
        e_SideChainVanDerWaalsSurfaceArea,        //!< Average van der waals surface area of the side chain
        e_SideChainHBondAcceptors,                //!< Number of hydrogen bond accepting atoms in the side chain
        e_SideChainHBondDonors,                   //!< Number of hydrogen bond donating atoms in the side chain (http://www.biochem.ucl.ac.uk/bsm/atlas/)
        e_Aromatic,                               //!< 1 for aromatic residues, 0 otherwise
        s_NumberPropertyTypes,                    //!< number of amino acid property types
        s_NumberPropertyTypesForANN = 7           //!< number of Properties that will be put into numerical representation of each AA for input to ANNs( jufo and contact)
      };

      //! @brief PropertyType as string
      //! @param PROPERTY_TYPE the PropertyType
      //! @return the string for the PropertyType
      static const std::string &GetPropertyDescriptor( const PropertyType &PROPERTY_TYPE);

      //! @brief PropertyTypeEnum is used for I/O of PropertyType
      typedef util::WrapperEnum< PropertyType, &GetPropertyDescriptor, s_NumberPropertyTypes> PropertyTypeEnum;

    private:

    //////////
    // data //
    //////////

      std::string                                                        m_ThreeLetterCode;                    //!< three letter code
      char                                                               m_OneLetterCode;                      //!< one letter code
      bool                                                               m_IsNaturalAminoAcid;                 //!< flag to indicate if the amino acid is one of the twenty natural occurring amino acids
      std::string                                                        m_ParentType;                         //!< aa type of parent amino acid for uncommon amino acids MSE would have MET
      storage::Set< AtomType>                                            m_AtomTypes;                          //!< atom types that are associated with this class
      storage::Map< AtomType, std::string>                               m_PDBAtomName;                        //!< atom name for within pdb files
      AtomType                                                           m_FirstSidechainAtomType;             //!< type of first side chain atom
      storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> > m_DihedralAtoms;                      //!< atoms that are used for calculating dihedral angles for the side chain of the amino acid. The atoms are ordered as needed to calculate dihedral angles starting at backbone and moving out along the side chain.
      AtomType                                                           m_SideChainIonizableHydrogen;         //!< Atom type that is ionizable at physiologically relevant pH (value given by other aa type data)
      double                                                             m_Properties[ s_NumberPropertyTypes]; //!< properties
      util::SiPtr< const AAType>                                         m_ParentTypePtr;                      //!< Parent type
      storage::Vector< chemistry::AtomType>                              m_ChemistryTypes;                     //!< Chemistry types
      storage::Vector< double>                                           m_VdwRadii;                           //!< Van-der waals radii for all biol::AtomTypes

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined amino acid
      AATypeData();

      //! @brief construct AATypeData with undefined properties
      //! @param THREE_LETTER_CODE Three letter code for this amino acid type
      //! @param ONE_LETTER_CODE One letter code for this amino acid type
      //! @param IS_NATURAL_AA boolean to indicate whether this amino acid type is one of 20 natural ones
      //! @param PARENT_TYPE name of parent AAType, this is the same as the AAType for natural amino acids
      AATypeData
      (
        const std::string &THREE_LETTER_CODE,
        const char ONE_LETTER_CODE,
        const bool IS_NATURAL_AA,
        const std::string &PARENT_TYPE
      );

      //! @brief construct amino acid from all its data
      AATypeData
      (
        const std::string &THREE_LETTER_CODE,
        const char ONE_LETTER_CODE,
        const bool IS_NATURAL_AA,
        const std::string &PARENT_TYPE,
        const storage::Set< AtomType> &ATOM_TYPES,
        const AtomType &FIRST_SIDECHAIN_ATOM_TYPE,
        const double NATURAL_PREVALENCE,
        const double STERICAL_PARAMETER,
        const double POLARIZABILITY,
        const double VOLUME,
        const double HYDROPHOBICITY,
        const double ISOELECTRIC_POINT,
        const double CHARGE,
        const double PK_EMBOSS,
        const double PK_DTASELECT,
        const double PK_SOLOMON,
        const double PK_SILLERO,
        const double PK_RODWELL,
        const double PK_PATRICKIOS,
        const double PK_WIKIPEDIA,
        const double PK_LEHNINGER,
        const double PK_GRIMSELY,
        const double PK_BJELLQVIST,
        const double PK_PROMOST,
        const double PK_BJELLQVIST_NTERM,
        const double PK_BJELLQVIST_CTERM,
        const double PK_PROMOST_NTERM,
        const double PK_PROMOST_CTERM,
        const double PK_CAREY_NTERM,
        const double PK_CAREY_CTERM,
        const double HELIX_PROBABILITY,
        const double STRAND_PROBABILITY,
        const double FREE_ENERGY_HELIX,
        const double FREE_ENERGY_STRAND,
        const double FREE_ENERGY_COIL,
        const double TRANSFER_FREE_ENERGY_WHIMLEY_WHITE,
        const double TRANSFER_FREE_ENERGY_ENGELMAN_SEITZ_GOLDMAN,
        const double TRANSFER_FREE_ENERGY_KYTE_DOOLITTLE,
        const double TRANSFER_FREE_ENERGY_EISENBERG,
        const double TRANSFER_FREE_ENERGY_HOPP_WOODS,
        const double TRANSFER_FREE_ENERGY_GUY,
        const double TRANSFER_FREE_ENERGY_JANIN,
        const double TRANSFER_FREE_ENERGY_PUNTA_MARITAN_1D,
        const double TRANSFER_FREE_ENERGY_PUNTA_MARITAN_3D,
        const double FREE_ENERGY_CORE,
        const double FREE_ENERGY_TRANSITION,
        const double FREE_ENERGY_SOLUTION,
        const double FREE_ENERGY_CORE_HELIX,
        const double FREE_ENERGY_TRANSITION_HELIX,
        const double FREE_ENERGY_SOLUTION_HELIX,
        const double FREE_ENERGY_CORE_STRAND,
        const double FREE_ENERGY_TRANSITION_STRAND,
        const double FREE_ENERGY_SOLUTION_STRAND,
        const double FREE_ENERGY_CORE_COIL,
        const double FREE_ENERGY_TRANSITION_COIL,
        const double FREE_ENERGY_SOLUTION_COIL,
        const double FREE_ENERGY_CORE_PORE,
        const double FREE_ENERGY_CORE_MEMBRANE,
        const double BETA_BARREL_HYDROPHOBICITY,
        const double FREE_ENERGY_EXT_BLAST_BB,
        const double FREE_ENERGY_EXT_TYPE_BB,
        const double SASA,
        const double SIDE_CHAIN_GIRTH,
        const double SIDE_CHAIN_POLARIZABILITY,
        const double SIDE_CHAIN_TOPOLOGICAL_POLAR_SURFACE_AREA,
        const double SIDE_CHAIN_VAN_DER_WAALS_SURFACE_AREA,
        const int    HBOND_ACCEPTORS,
        const int    HBOND_DONORS,
        const bool   IS_AROMATIC,
        const double VDW_RADIUS_CB
      );

      //! @brief virtual copy constructor
      AATypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get three letter code
      //! @return three letter code
      const std::string &GetThreeLetterCode() const
      {
        return m_ThreeLetterCode;
      }

      //! @brief get one letter code
      //! @return one letter code
      char GetOneLetterCode() const
      {
        return m_OneLetterCode;
      }

      //! @brief get the parent amino acid type
      //! this maps uncommon aas to their natural counterpart like Seleno-methionin to methionin
      //! @return AAType of parent for this uncommon amino acid
      const AAType &GetParentType() const;

      //! @brief get the parent amino acid type threelettercode
      //! this maps uncommon aas to their natural counterpart like Seleno-methionin to methionin
      //! @return three letter code of AAType of parent for this uncommon amino acid
      const std::string &GetParentTypeString() const
      {
        return m_ParentType;
      }

      //! @brief get if amino acid is a natural amino acid
      //! @return true, if this is data for one of the twnety natural amino acids
      bool IsNaturalAminoAcid() const
      {
        return m_IsNaturalAminoAcid;
      }

      //! @brief atom types for that amino acid type
      //! @return a set of atom types that are part of that amino acid type
      const storage::Set< AtomType> &GetAllowedAtomTypes() const
      {
        return m_AtomTypes;
      }

      //! @brief atom types for that amino acid type without hydrogens
      //! @return a set of atom types that are part of that amino acid type without hydrogens
      storage::Set< AtomType> GetAllowedHeavyAtomTypes() const;

      //! @brief return the pdb atom name for the given AtomType
      //! @param ATOM_TYPE atom type
      //! @return the pdb atom name with proper spacing, so that it can be written to the pdb file, empty if this amino acid type does not have this atom type
      const std::string &GetPDBAtomName( const AtomType &ATOM_TYPE) const;

      //! @brief return the associated atom type from a given atom name (IUPAC or PDB)
      //! @param ATOM_NAME the atom name defined by either the pdb or IUPAC for this aa type
      //! @return the atom type for that pdb atom name, if it is defined for that amino acid type
      AtomType GetAtomTypeFromAtomName( const std::string &ATOM_NAME) const;

      //! @brief type of first side chain atom
      //! @return type of first side chain atom - HA2 for GLY or similar, CB for most other Amino acids
      const AtomType &GetFirstSidechainAtomType() const
      {
        return m_FirstSidechainAtomType;
      }

      //! @brief get requested PROPERTY
      //! @param PROPERTY AAProperty of interet
      //! @return PROPERTY
      double GetAAProperty( const PropertyType PROPERTY) const
      {
        return m_Properties[ PROPERTY];
      }

      //! @brief atoms that are used for calculating dihedral angles for the side chain of the amino acid.
      //! The atoms are ordered as needed to calculate dihedral angles starting at backbone and moving out
      //! along the side chain.
      //! @return storage vector of atom types
      storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> > GetSideChainDihedralAngleAtomTypes() const
      {
        return m_DihedralAtoms;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief map a given atom type to a pdb atom name - is only used in the constructor of the AtomTypes class
      //! @param ATOM_TYPE the atom type
      //! @param PDB_ATOM_NAME the atom name defined by the pdb for this aa type
      //! @return true if mapping was successful (false if the atom type was already added)
      bool AddAtomTypePDBAtomNameMapping( const AtomType &ATOM_TYPE, const std::string &PDB_ATOM_NAME);

      //! @brief check if given ATOM_TYPE is part of the amino acid
      //! @param ATOM_TYPE type that is checked for
      //! @return true is this amino acid type has the given ATOM_TYPE
      bool DoesContainAtomType( const AtomType &ATOM_TYPE) const;

      //! @brief get properties to be used in ANN as a vector
      //! @return properties to be used in ANN as a vector
      linal::Vector< double> GetPropertiesForANN() const;

      //! @brief get a chemistry::FragmentConfigurationalShared representation of the AA
      //! @param C_TERMINAL whether to return a representation with the C-terminal carboxyl group
      //! @param N_TERMINAL whether to return a representation with the N-terminal amino group
      //! @return a chemistry::FragmentConfigurationalShared representation of the AA
      const chemistry::FragmentConfigurationShared &GetFragment
      (
        const bool &C_TERMINAL,
        const bool &N_TERMINAL
      ) const;

      //! @brief Get the chemistry atom type of a particular atom in this aa type
      //! @param ATOM the biol::AtomType of interest
      //! @return chemistry::AtomType of that Atom
      chemistry::AtomType GetChemistryAtomType( const AtomType &ATOM_TYPE) const;

      //! @brief Get the effective VdW radius of an atom type for this AA to other AAs
      //! @param ATOM the biol::AtomType of interest
      //! @return van-der-waals radius to external AAs of the given atom type
      double GetVdwRadiusToOtherAA( const AtomType &ATOM_TYPE) const;

      //! @brief gets the structure factor for this AA Type
      //! @return the structure factor for this AA Type
      util::ShPtr< math::FunctionInterfaceSerializable< restraint::SasDataParameters, double> > GetStructureFactor() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief set the dihedral angles for this AA type
      //! @param ANGLES the dihedral angles for this AAType
      void SetDihedralAngles( const storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> > &ANGLES);

      //! @brief set the ionizable H
      //! @param ATOM_TYPE the ionizable hydrogen
      void SetSideChainIonizableHType( const AtomType &ATOM_TYPE);

      //! @brief get the counts of what type and how many atoms in the side chain of a given amino acid
      //! @return Return the map of type and count of atoms in the side chain of a given amino acid
      storage::Map< chemistry::ElementType, size_t> GetElementTypeCount() const;

      //! @brief Calculate the mass of this amino acid (not when peptide bonded!)
      //! @return the mass of the amino acid as the sum of all atom masses; non-monoisotopic
      double CalculateMass() const;

      //! @brief determine the parent type of amino acid
      void DetermineParentType( const AATypes &AATYPES);

    }; //class AATypeData

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_TYPE_DATA_H_
