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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_ISOLATE_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_ISOLATE_H_

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
    //! @class FragmentSplitIsolate
    //! @brief This class get all components of a molecule complex that is provided.
    //!
    //! @see @link example_chemistry_fragment_split_isolate.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Oct 23, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitIsolate :
      public FragmentSplitInterface
    {

    private:

        size_t      m_MinSize;    //!< minimum size of the components that are desired

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      FragmentSplitIsolate *Clone() const;

      //! brief constructor
      FragmentSplitIsolate( const size_t MIN_SIZE = size_t( 2));

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
      const size_t GetMinSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns connected components of a graph that is not connected
      //! @param MOLECULE molecule of interest
      //! @param MOLECULE_GRAPH graph of molecule of interest
      //! @return connected components of a graph that is not connected
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

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_ISOLATE_H_
