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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_LINEAR_FRAGMENTS_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_LINEAR_FRAGMENTS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_split_interface.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentSplitLinearFragments
    //! @brief This class splits a molecule into fragments by growing them out one atom at a time in a linear fashion
    //! @brief similar to the openbabel FP2 fingerprint (http://openbabel.org/wiki/FP2)
    //!
    //! @see @link example_chemistry_fragment_split_linear_fragments.cpp @endlink
    //! @author geanesar
    //! @date Sep 22, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitLinearFragments :
      public FragmentSplitInterface
    {

    private:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! the maximum number of bonds to extend for each fragment
      size_t m_Steps;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      FragmentSplitLinearFragments( const size_t &STEPS = 7);

      //! @brief copy constructor
      //! @return a pointer to a copy of this class
      FragmentSplitLinearFragments *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief Get a description for what this class does (used when writing help)
      //! @return a description for what this class does (used when writing help)
      const std::string &GetClassDescription() const;

      //! @brief gets the minimum size of fragments
      //! @return the minimum size of fragments
      const size_t GetMinSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief a recursive algorithm for getting the fragments of each molecule
      //! @param MOLECULE_GRAPH the atom graph of the molecule of interest
      //! @param CURRENT_CHAIN a set of indices that represent the current linear fragment up to the current vertex
      //! @param CURRENT_INDEX the index that should be added to CURRENT_CHAIN, if possible
      //! @param FRAGMENTS a set of vectors that represent all unique linear fragments that have been seen
      //! @param MAX_STEPS the maximum number of bonds to go out
      //! @param IGNORE_H whether to ignore hydrogens or not
      void GetFragments
      (
        const ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH,
        const storage::Set< size_t> &CURRENT_CHAIN,
        const size_t &CURRENT_INDEX,
        storage::Set< storage::Vector< size_t> > &FRAGMENTS,
        const size_t &MAX_STEPS,
        const bool &IGNORE_H = false
      ) const;

      //! @brief returns list of ring or chain fragments
      //! @param MOLECULE molecule of interest
      //! @param MOLECULE_GRAPH graph of molecule of interest
      //! @return list of ring or chain fragments
      storage::List< storage::Vector< size_t> > GetComponentVertices
      (
        const ConformationInterface &MOLECULE,
        ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_LINEAR_FRAGMENTS_H_
