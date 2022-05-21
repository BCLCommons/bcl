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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_ECFP_FRAGMENTS_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_ECFP_FRAGMENTS_H_

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
    //! @class FragmentSplitECFPFragments
    //! @brief This class splits a molecule into fragments by growing them out one atom at a time, similar to how
    //! @brief ECFP fragments are generated
    //!
    //! @details for additional information @see @link http://pubs.acs.org/doi/abs/10.1021/ci100050t @endlink
    //! @see @link example_chemistry_fragment_split_ecfp_fragments.cpp @endlink
    //! @author geanesar
    //! @date Sep 19, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitECFPFragments :
      public FragmentSplitInterface
    {

    private:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! the number of bonds from an atom to extend each fragment
      size_t m_Steps;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      FragmentSplitECFPFragments( const size_t &STEPS = 4);

      //! @brief copy constructor
      //! @return a pointer to a copy of this class
      FragmentSplitECFPFragments *Clone() const;

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

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_ECFP_FRAGMENTS_H_
