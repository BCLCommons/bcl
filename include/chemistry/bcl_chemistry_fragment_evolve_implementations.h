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

#ifndef BCL_CHEMISTRY_FRAGMENT_EVOLVE_IMPLEMENTATIONS_H_
#define BCL_CHEMISTRY_FRAGMENT_EVOLVE_IMPLEMENTATIONS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_evolve_base.h"
#include "bcl_chemistry_reaction_ensemble.h"
#include "bcl_chemistry_reaction_worker.h"
#include "io/bcl_io_serialization.h"
#include "sdf/bcl_sdf_rxn_factory.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentEvolveImplementations
    //! @brief Class used for mutating FragmentCompletes by mixing two FragmentCompletes together or mutating them
    //!
    //! @see @link example_chemistry_fragment_evolve_implementations.cpp @endlink
    //! @author geanesar, brownbp1
    //! @date Jan 14, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API FragmentEvolveImplementations :
      public util::SerializableInterface
    {
    public:

      //! Enum for specifying what kind of evolution to use
      enum EvolveType
      {
        e_Clone = 0,
        e_FragAdd,
        e_FragDel,
        e_Combine,
        s_End
      };

    private:

    //////////
    // data //
    //////////

      //! The type of evolution operation this class will do
      EvolveType m_EvolveType;

      //! The filename to read fragments from (if applicable)
      std::string m_FragmentFilename;

      //! The fragment pool
      mutable util::ShPtr< FragmentEnsemble> m_FragmentPool;
      
      //! Whether to randomize the return if multiple molecules are generated
      bool m_Randomize;

      //! The maximum size of concerned fragments
      size_t m_MaxNumberAtoms;

      //! Whether to check if bad bonds were formed
      bool m_AllowBadBonds;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor specifying evolve type
      //! @param EVOLVE_TYPE the evolution type to use
      //! @param CHECK_BAD_BONDS whether to ensure no bad bonds were formed
      FragmentEvolveImplementations( const EvolveType &EVOLVE_TYPE, const bool &ALLOW_BAD_BONDS = true);

      //! @brief Clone function
      //! @return pointer to new FragmentEvolveImplementations
      FragmentEvolveImplementations *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns class name of the object when used in a dynamic context
      //! @return the class name
      const std::string &GetAlias() const;

      //! @brief gets the type of evolution operations that will be performed
      //! @return the evolution type
      const EvolveType &GetEvolveType() const;

      //! @brief gets the string associated with each evolve type
      //! @param EVOLVE_TYPE the evolve type to get the string for
      //! @return a string describing the evolve type
      static const std::string &GetEvolveTypeString( const FragmentEvolveImplementations::EvolveType &EVOLVE_TYPE);

      //! @brief gets the fragment list that is in use
      //! @return a ShPtr to a FragmentEnsemble that holds the fragments
      const util::ShPtr< FragmentEnsemble> &GetFragmentList() const;

      //! @brief sets the fragment list for this class
      //! @param FRAGMENT_LIST the fragment list to use
      void SetFragmentList( const FragmentEnsemble &FRAGMENT_LIST);

      //! @brief sets the fragment list for this class
      //! @param FRAGMENT_LIST the fragment list to use
      void SetFragmentList( const util::ShPtr< FragmentEnsemble> &FRAGMENT_LIST);

      //! @brief reads in fragments from a file and sets up the fragment list for this class
      //! @param FILENAME the file to read from
      void SetFragmentList( const std::string &FILENAME);

      //! @brief sets whether to check for bad bonds
      //! @param ALLOW_BAD_BONDS if false, adding fragments will check for bad bonds
      void SetAllowBadBonds( const bool &ALLOW_BAD_BONDS);

      //! @brief gets the number of molecules needed by operator()
      //! @return the number needed
      size_t NumRequiredMols() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief run the operation specified by the class
      //! @param MEMBERS an ensemble of molecules to use
      //! @return a new ensemble of mutated molecules
      util::ShPtrVector< FragmentComplete> MakeMolecules
      (
        const util::SiPtrVector< const FragmentComplete> &MEMBERS
      ) const;

      //! @brief combines sections of two molecules
      //! @param FIRST_MOLECULE the first molecule
      //! @param SECOND_MOLECULE the second molecule
      //! @param MAX_FRAG_SIZE the maximum size of fragments that should be connected together
      //! @return a new fragment consisting of pieces of the inputs
      util::ShPtrVector< FragmentComplete> Combine
      (
        const FragmentComplete &FIRST_MOLECULE,
        const FragmentComplete &SECOND_MOLECULE,
        const size_t MAX_FRAG_SIZE = 50
      ) const;

      //! @brief add a fragment to a molecule
      //! @param MOLECULE the molecule to add a fragment to
      //! @return a new molecule with an added fragment
      util::ShPtrVector< FragmentComplete> MutateAdd
      (
        const FragmentComplete &MOLECULE,
        const FragmentEnsemble &FRAGMENTS
      ) const;

      //! @brief remove a fragment from a molecule
      //! @param MOLECULE the molecule to remove a fragment from
      //! @return fragments of the molecule resulting from
      util::ShPtrVector< FragmentComplete> MutateDel( const FragmentComplete &MOLECULE) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief determines what fragments would result from breaking a bond in a graph
      //! @param MOLECULE_GRAPH the graph that will have its bond broken
      //! @param FROM one vertex
      //! @param TO the other vertex
      //! @return a list of vectors of indices which correspond to connected components of the graph
      storage::List< storage::Vector< size_t> > CollectFragmentsFromBondBreakage
      (
        graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
        const size_t &FROM,
        const size_t &TO
      ) const;

      //! @brief determines what fragments would result from breaking a bond in a graph
      //! @param MOLECULE the molecule that will be fragmented
      //! @param MOLECULE_GRAPH the graph that will have its bond broken
      storage::List< storage::Vector< size_t> > FragmentsFromRandomBondBreakage
      (
        const FragmentComplete &MOLECULE,
        graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
        const size_t &EDGE_TYPE = 1
      ) const;

    ////////////////
    // operations //
    ////////////////
    
      //! @brief reads in fragments from a file
      //! @param FILENAME the file to read fragments from
      //! @return true if the fragments were read
      bool ReadFragmentsFromFile() const;

      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &SERIALIZER,
        std::ostream &ERR_STREAM
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class FragmentEvolveImplementations

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_EVOLVE_IMPLEMENTATIONS_H_
