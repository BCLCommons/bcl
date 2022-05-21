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
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "fold/bcl_fold_placement_domain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> PlacementDomain::s_Instance
    (
      GetObjectInstances().AddInstance( new PlacementDomain())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PlacementDomain::PlacementDomain()
    {
    }

    //! @brief returns a pointer to a new PlacementDomain
    //! @return pointer to a new PlacementDomain
    PlacementDomain *PlacementDomain::Clone() const
    {
      return new PlacementDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &PlacementDomain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief computes a transformation matrix for the placement of the given domain in the given protein model
    //! @param SELECTED_DOMAIN domain which shall be placed in the protein model
    //! @param PROTEIN_MODEL protein model in which the domain shall be placed
    //! @return the matrix defining the transformation to be applied on the given domain to place it in the given
    //! protein model and a boolean value indicating if the computation of the transformation matrix was successful
    storage::Pair< math::TransformationMatrix3D, bool> PlacementDomain::Place
    (
      const assemble::Domain &SELECTED_DOMAIN,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // holds the SSE in which neighborhood to place the domain
      util::SiPtr< const assemble::SSE> p_neighbor_sse;

      // find the SSE in the model with the smallest sequence distance to the domain TODO
      util::SiPtrVector< const assemble::SSE> tmp( PROTEIN_MODEL.GetSSEs());
      util::SiPtrList< const assemble::SSE> sses_in_model( tmp.Begin(), tmp.End());

      // if no suitable SSE could be found pick a random SSE
      if( !p_neighbor_sse.IsDefined())
      {
        p_neighbor_sse = assemble::PickSSERandom().Pick( sses_in_model);
      }

      // compute the transformation matrix for the starting from the orientation of the selected SSE
      math::TransformationMatrix3D transformation;
      if( p_neighbor_sse.IsDefined())
      {
        transformation( p_neighbor_sse->GetOrientation());
      }

      // add random transformations TODO add randomized flipping etc
      const random::DistributionInterface &random( random::GetGlobalRandom());
      transformation
      (
        random.Random< double>( 0.0, 10.0), random.Random< double>( 0.0, 10.0), random.Random< double>( 0.0, 10.0)
      );

      return storage::Pair< math::TransformationMatrix3D, bool>( transformation, true);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read object from input stream
    //! @param ISTREAM input stream to read object from
    //! @return input stream which was read from
    std::istream &PlacementDomain::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write object into  output stream
    //! @param OSTREAM output stream to write object into
    //! @param INDENT number of indentations to separate members
    //! @return output stream object was written into
    std::ostream &PlacementDomain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
