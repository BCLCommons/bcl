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
#include "fold/bcl_fold_mutate_domain_transformation.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDomainTransformation::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDomainTransformation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDomainTransformation::MutateDomainTransformation() :
      m_Transformer()
    {
    }

    //! @brief constructor taking parameters
    //! @param TRANSFORMER method to transform the transformation matrix 3d of the domain
    MutateDomainTransformation::MutateDomainTransformation
    (
      const util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > &TRANSFORMER
    ) :
      m_Transformer( TRANSFORMER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDomainTransformation
    MutateDomainTransformation *MutateDomainTransformation::Clone() const
    {
      return new MutateDomainTransformation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDomainTransformation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param MUTATE_DOMAIN the domain that is going to be mutated
    //! @return mutate result pointing to mutated domain
    math::MutateResult< assemble::Domain>
    MutateDomainTransformation::operator()( const assemble::Domain &MUTATE_DOMAIN) const
    {
      // get current orientation of the domain
      const math::TransformationMatrix3D current_orientation( MUTATE_DOMAIN.GetOrientation());

      // mutate the current orientation
      const math::MutateResult< math::TransformationMatrix3D> mutated_orientation
      (
        m_Transformer->operator()( current_orientation)
      );

      // make new domain copied from old domain
      util::ShPtr< assemble::Domain> new_domain( MUTATE_DOMAIN.Clone());

      // set the new domain to the the new orientation
      new_domain->Transform( *mutated_orientation.GetArgument());

      // make mutate result out of mutated domain
      const math::MutateResult< assemble::Domain> mutate_result( new_domain, *this);

      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDomainTransformation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Transformer, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDomainTransformation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Transformer, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
  
} // namespace bcl
