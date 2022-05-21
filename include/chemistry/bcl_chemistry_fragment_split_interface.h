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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_INTERFACE_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "util/bcl_util_function_interface_serializable.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentSplitInterface
    //! @brief This class is an interface class for classes that split molecules
    //!
    //! @see @link example_chemistry_fragment_split_interface.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Oct 23, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitInterface :
      public virtual util::FunctionInterfaceSerializable< ConformationInterface, FragmentEnsemble>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      virtual FragmentSplitInterface *Clone() const = 0;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      virtual const std::string &GetAlias() const = 0;

      //! get the minimum size of a component of interest
      virtual const size_t GetMinSize() const = 0;

      //! @brief returns connected components of a graph that is not connected
      //! @param MOLECULE molecule of interest
      //! @param MOLECULE_GRAPH graph of molecule of interest
      //! @return connected components of a graph that is not connected
      virtual storage::List< storage::Vector< size_t> > GetComponentVertices
      (
        const ConformationInterface &MOLECULE,
        ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
      ) const
      {
        BCL_Exit( "Class has not implemented GetComponentVertices; use operator() instead!", -1);
        return storage::List< storage::Vector< size_t> >();
      }

      //! @brief splits the molecule according to GetComponentVertices
      //! @param CONFORMATION the molecule to split
      virtual FragmentEnsemble operator()( const ConformationInterface &CONFORMATION) const;

      //! @brief helper function that can be used to convert the output of GetComponentVertices into an ensemble
      //! @param MOLECULE the molecule to split
      //! @param COMPONENTS list of components to create
      //! @param MOLECULE_GRAPH graph of the molecule with atoms
      //! @return a fragment ensemble
      virtual FragmentEnsemble ConvertComponentsIntoEnsemble
      (
        const ConformationInterface &MOLECULE,
        const storage::List< storage::Vector< size_t> > &COMPONENTS,
        const ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH,
        const bool &INVERT_OUTPUT = false
      ) const;

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_INTERFACE_H_
