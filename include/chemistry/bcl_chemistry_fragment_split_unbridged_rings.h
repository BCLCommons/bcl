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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_UNBRIDGED_RINGS_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_UNBRIDGED_RINGS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_split_interface.h"
#include "util/bcl_util_enumerated.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentSplitUnbridgedRings
    //! @brief This class returns rings or chains for a molecule that is provided
    //!
    //! @see @link example_chemistry_fragment_split_unbridged_rings.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Oct 23, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitUnbridgedRings :
      public FragmentSplitInterface
    {

    private:

        bool      m_Aromatic;

        size_t    m_MaxSize;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_InstanceRing;
      static const util::SiPtr< const util::ObjectInterface> s_InstanceChain;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      FragmentSplitUnbridgedRings *Clone() const;

      //! @brief constructor
      //! @param GET_RINGS true if rings are desired, false if chains are desired
      //! @param MIN_SIZE get the minimum size of ring or chain that is desired
      FragmentSplitUnbridgedRings( bool GET_AROMATIC, const size_t MAX_SIZE = size_t( 100));

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

      //! get the minimum size of a component of interest
      const size_t GetMaxSize() const;

      //! get the minimum size of a component of interest
      const size_t GetMinSize() const
      {
        return size_t( 2);
      }

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

      static bool IsBridgeHead( const AtomConformationalInterface &ATOM);

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

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_UNBRIDGED_RINGS_H_
