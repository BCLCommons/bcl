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

#ifndef BCL_FOLD_MUTATE_DOMAIN_TRANSFORMATION_H_
#define BCL_FOLD_MUTATE_DOMAIN_TRANSFORMATION_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDomainTransformation
    //! @brief mutates the orientation of a domain according to a transformation matrix 3d mutate object
    //! @details uses an object that mutates transformation matrix 3d objects to change the orientation of a domain
    //!
    //! @see @link example_fold_mutate_domain_transformation.cpp @endlink
    //! @author alexanns
    //! @date Apr 20, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDomainTransformation :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! method to transform the transformation matrix 3d of the domain
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > m_Transformer;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateDomainTransformation();

      //! @brief constructor taking parameters
      //! @param TRANSFORMER method to transform the transformation matrix 3d of the domain
      MutateDomainTransformation( const util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > &TRANSFORMER);

      //! @brief Clone function
      //! @return pointer to new MutateDomainRotationTranslation
      MutateDomainTransformation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param MUTATE_DOMAIN the domain that is going to be mutated
      //! @return mutate result pointing to mutated domain
      math::MutateResult< assemble::Domain> operator()( const assemble::Domain &MUTATE_DOMAIN) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class MutateDomainTransformation

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_DOMAIN_TRANSFORMATION_H_ 
