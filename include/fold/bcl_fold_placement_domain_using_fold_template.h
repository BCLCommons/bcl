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

#ifndef BCL_FOLD_PLACEMENT_DOMAIN_USING_FOLD_TEMPLATE_H_
#define BCL_FOLD_PLACEMENT_DOMAIN_USING_FOLD_TEMPLATE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_domain_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementDomainUsingFoldTemplate
    //! @brief places SSEs in domain into a protein model using a fold template
    //!
    //! @see @link example_fold_placement_domain_using_fold_template.cpp @endlink
    //! @author weinerbe
    //! @date Feb 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementDomainUsingFoldTemplate :
      public PlacementDomainInterface
    {

    private:

    //////////
    // data //
    //////////

      //! body deviation; positive values get a larger template, and negative values exclude SSEs from the pool
      int m_BodyDeviation;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor with body deviations
      //! @param BODY_DEVIATION positive values get a larger template, and negative values exclude SSEs
      PlacementDomainUsingFoldTemplate( const int &BODY_DEVIATION = 0);

      //! @brief Clone function
      //! @return pointer to new PlacementDomainUsingFoldTemplate
      PlacementDomainUsingFoldTemplate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief determines transformation matrices for placing SSEs from a Domain into a protein model
      //! @param DOMAIN_TO_PLACE domain containing SSEs to be placed
      //! @param PROTEIN_MODEL model that domain will be placed into
      //! @return transformation matrices for placing SSEs from a Domain into a protein model
      storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool> Place
      (
        const assemble::DomainInterface &DOMAIN_TO_PLACE, const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief determines transformation matrices for placing SSEs from a Domain into a protein model
      //! @param DOMAIN_TO_PLACE domain containing SSEs to be placed
      //! @return transformation matrices for placing SSEs from a Domain into a protein model
      storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool> Place
      (
        const assemble::DomainInterface &DOMAIN_TO_PLACE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PlacementDomainUsingFoldTemplate

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PLACEMENT_DOMAIN_USING_FOLD_TEMPLATE_H_ 
