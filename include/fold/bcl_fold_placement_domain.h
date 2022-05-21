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

#ifndef BCL_FOLD_PLACEMENT_DOMAIN_H_
#define BCL_FOLD_PLACEMENT_DOMAIN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementDomain
    //! @brief computes a transformation matrix to place a domain in a protein model
    //! @detail this class tries to place the domain in the neighborhood of a SSE in the protein model
    //!
    //! @author fischea
    //! @date April 19, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementDomain :
      public PlacementInterface< assemble::Domain, assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PlacementDomain();

      //! @brief returns a pointer to a new PlacementDomain
      //! @return pointer to a new PlacementDomain
      PlacementDomain *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief computes a transformation matrix for the placement of the given domain in the given protein model
      //! @param SELECTED_DOMAIN domain which shall be placed in the protein model
      //! @param PROTEIN_MODEL protein model in which the domain shall be placed
      //! @return the matrix defining the transformation to be applied on the given domain to place it in the given
      //! protein model and a boolean value indicating if the computation of the transformation matrix was successful
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::Domain &SELECTED_DOMAIN,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read object from input stream
      //! @param ISTREAM input stream to read object from
      //! @return input stream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write object into  output stream
      //! @param OSTREAM output stream to write object into
      //! @param INDENT number of indentations to separate members
      //! @return output stream object was written into
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PlacementDomain

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PLACEMENT_DOMAIN_H_
