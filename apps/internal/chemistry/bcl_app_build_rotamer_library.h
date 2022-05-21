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
#ifndef BCL_APP_BUILD_ROTAMER_LIBRARY_H_
#define BCL_APP_BUILD_ROTAMER_LIBRARY_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "descriptor/bcl_descriptor.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// include header of this class
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BuildRotamerLibrary
    //! @brief Application for creating rotamer library for fragments using conformations from a structure database
    //!
    //! @see @link example_app_build_rotamer_library.cpp @endlink
    //! @author kothiwsk
    //! @date 10/29/2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BuildRotamerLibrary :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! database/file containing molecule structures from which to find fragments
      util::ShPtr< command::FlagInterface>                                          m_StructureDatabaseFlag;

      //! flag for specifying if we want fragments so as to cover all dihedral bond types
      util::ShPtr< command::FlagInterface>                                          m_DihedralFragmentFlag;

      //! flag for specifying if we want fragments so as to cover all bond angle types
      util::ShPtr< command::FlagInterface>                                          m_BondAngleFragmentFlag;

      //! flag for specifying if rotamers need to be updated in a given rotamer libray, given a structure database
      util::ShPtr< command::FlagInterface>                                          m_IdentifyNewConformersFlag;

      //! the fragment will only be considered a novel conformation if the conformations returns a # below this when
      //! run through the comparer
      util::ShPtr< command::FlagInterface>                                          m_BinSizeFlag;

      //! output filename to store rotamer library
      util::ShPtr< command::FlagInterface>                                          m_OutputFileFlag;

      //! minimum fragment counts required for saving rotamer library of a particular fragment
      util::ShPtr< command::FlagInterface>                                          m_MinFragmentCount;

      //! maximum fragment counts for fragment to be considered for finding cluster center. Larger the counts more the
      //! time it takes to find cluster center
      util::ShPtr< command::FlagInterface>                                          m_MaxFragmentCount;

      //! minimum rotamer counts required for saving particular rotamer of a fragment
      util::ShPtr< command::FlagInterface>                                          m_MinRotamerCount;

      //! weight to be given to new rotamers seen in the conformer library
      util::ShPtr< command::FlagInterface>                                          m_CountFactor;

      //! weight to be given to new rotamers seen in the conformer library
      util::ShPtr< command::FlagInterface>                                          m_GetIndividualRotamers;

      //! option to add simulated counts for all chain dihedrals
      util::ShPtr< command::FlagInterface>                                          m_AddSimulatedCounts;

      //! object that matches whole ring like searching for whole word against an incomplete word
      mutable util::ShPtr< chemistry::FragmentGraphMarker>                          m_MapRings;

      //! simple graphs of all the molecules in structure database
      mutable storage::Vector< graph::ConstGraph< size_t, size_t> >                 m_MoleculeSimpleGraphs;

      //! simple graphs of all the fragments
      mutable storage::Vector< graph::ConstGraph< size_t, size_t> >                 m_FragmentSimpleGraphs;

      //! atom graphs of all the molecules in the structure databse
      mutable storage::Vector
      <
        graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t>
      >                                                                             m_EnsembleGraphs;

      //! the output file, which writes out the rotamer library of fragments as they are generated
      mutable io::OFStream                                                          m_Output;

      //! bin size to be used for comparing conformations
      mutable double                                                                m_BinSize;

      //! Map to hold unique dihedral angle types in string format and the molecule
      mutable storage::Map< std::string, chemistry::FragmentComplete>               m_DihedralMap;

      // Key <- atom-type <> Vector< bond type, atom type>
      // Value <- Matrix with unit-vector coordinates of all atoms after the first.
      // The first atom is always moved to 1 0 0, second atom is moved such that it is 0 in
      mutable storage::Map
      <
        typename chemistry::RotamerLibraryInterface::t_BondAngleMapKey,
        storage::List< storage::Pair< linal::Matrix< double>, math::RunningAverage< linal::Vector< double> > > >
      > m_BondAngles;

      //! conformation sampler, if m_AddSimulatedCounts is turned on
      mutable chemistry::SampleConformations m_SampleConformations;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      BuildRotamerLibrary();

      //! copy constructor
      BuildRotamerLibrary( const BuildRotamerLibrary &APP);

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType BuildRotamerLibrary_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      BuildRotamerLibrary *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the full description
      //! @return the full description
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetWebText() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief determine rings and bonds contained in the fragment
      void FindDihedralFragments() const;

      //! @brief determine rings and bonds contained in the fragment
      void FindBondAngleFragments() const;

      //! @brief determine rings and bonds contained in the fragment
      //! @param FRAGMENT molecule whose rings are to be determined
      //! @return a list of list of indices of atoms. Inner list has atom indices contained in a single ring.
      storage::List< storage::List< size_t> > DetermineRings( const chemistry::FragmentComplete &FRAGMENT) const;

      //! @brief determine if rings of molecule contain multiple ring rotamers
      //! @param RING_BONDS a list of list of indices of atoms. Inner list has atom indices contained in a single ring
      //! @param BIN_VECTOR rotamer bins for the corresponding rings
      //! @return true if any ring has multiple rotamers else false
      bool IfContainsRingConformations
      (
        const storage::List< storage::List< size_t> > &RING_BONDS,
        const storage::List< storage::Vector< size_t> > &BIN_VECTOR
      ) const;

      //! @brief finds all the conformers of a fragment inside an ensemble
      //! @param FRAGMENT the fragment to find conformers of
      //! @param FRAGMENT_INDEX index of the fragment in the fragment ensemble
      //! @return void
      void FindConformers
      (
        const chemistry::FragmentComplete &FRAGMENT,
        const size_t &FRAGMENT_INDEX
      ) const;

      //! @brief updates rotamer information for given fragment given a set of rotamers in the form of RotamerClusterCenter
      //! @param FRAGMENT the fragment whose rotamers need to be updated
      //! @param FRAGMENT_INDEX index of the fragment in the fragment ensemble
      //! @param CLUSTER_CENTER contains rotamer information that needs to be added to the fragment of interest
      //! @param RING_BONDS ring bonds of fragment.
      //! @return void
      void AddRotamers
      (
        const chemistry::FragmentComplete &FRAGMENT,
        const size_t FRAGMENT_INDEX,
        const chemistry::RotamerClusterCenter &CLUSTER_CENTER,
        const storage::List< storage::List< size_t> > &RING_BONDS
      ) const;

      //! @brief determine wheter rotamers have different ring conformations
      //! @param FRAGMENT rotamers for a fragment of interest
      //! @param FRAGMENT_INDEX rotamers that have been determined for the fragment of interest
      //! @return true if different ring conformations exist, otherwise false
      bool HasRingConformations
      (
        const chemistry::RotamerEnsemble &FRAGMENT,
        const storage::List< linal::Vector< int> > &FRAGMENT_ROTAMERS
      ) const;

      //! @brief finds all the conformers of a fragment inside an ensemble
      //! @param FRAGMENT the fragment to find conformers of
      //! @param FRAGMENT_INDEX index of the fragment in the fragment ensemble
      //! @param CLUSTER_CENTER object containing all rotamers for the fragment of interest
      //! @param RING_BONDS true if fragment contains ring, false if not.
      void OutputRotamerLibrary
      (
        const chemistry::FragmentComplete &FRAGMENT,
        const size_t FRAGMENT_INDEX,
        const chemistry::RotamerClusterCenter &CLUSTER_CENTER,
        const storage::List< storage::List< size_t> > &RING_BONDS
      ) const;

      //! @brief convert dihedral bond information into a string
      //! @param TYPE dihedral bond information
      template< unsigned int t_N>
      std::string ConvertAtomBondTypeIntoString( const storage::VectorND< t_N, size_t> &TYPE) const;

      //! @brief convert dihedral information into molecule
      //! @TYPE dihedral bond information
      //! @return molecule from dihedral bond
      chemistry::FragmentComplete ConvertDihedralToMolecule( const storage::VectorND< 7, size_t> &TYPE) const;

      //! @brief get dihedral bonds as fragments from the given molecule
      //! @param MOLECULE the small molecule from which dihedral information will be obtained
      void GetDihedralFragments( const chemistry::ConformationInterface &MOLECULE) const;

      //! @brief get dihedral bonds as fragments from the given molecule
      //! @param MOLECULE the small molecule from which dihedral information will be obtained
      void GetBondAngleFragments( const chemistry::FragmentComplete &MOLECULE) const;

    }; // BuildRotamerLibrary

  } // namespace app
} // namespace bcl
#endif // BCL_APP_BUILD_ROTAMER_LIBRARY_H_
