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

#ifndef BCL_FOLD_MUTATE_DOMAIN_FLIP_H_
#define BCL_FOLD_MUTATE_DOMAIN_FLIP_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDomainFlip
    //! @brief flips a domain  externally or internally (each SSE by itself)
    //! @details This Mutate apply a flip, 180 degree rotation around a given axis, internally or externally. Another boolean
    //! decides whether to flip all strands or a randomly selected subset. The external rotation should only be used
    //! when coupled with SheetCollectors, since the Domain itself does not defined a valid rotation.
    //!
    //! @see @link example_fold_mutate_domain_flip.cpp @endlink
    //! @author karakam
    //! @date Mar 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDomainFlip :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! flip axes
      storage::Set< coord::Axis> m_FlipAxes;

      //! boolean to flip internally
      bool m_FlipInternal;

      //! boolean to decide whether to flip all strands or only few
      bool m_FlipAll;

      //! boolean to decide whether to use different flip axes for each flip, valid if flip internal
      bool m_UseDifferentFlipAxes;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a single flip axis, flip internal and and a flip all boolean
      //! @param FLIP_AXIS Axis to apply the flip around
      //! @param FLIP_INTERNAL boolean to flip internally
      //! @param FLIP_ALL boolean to decide whether to flip all strands or only few
      //! @param USE_DIFFERENT_FLIP_AXES boolean to decide whether to use different flip axes for each flip
      MutateDomainFlip
      (
        const coord::Axis FLIP_AXIS = coord::GetAxes().e_Z,
        const bool FLIP_INTERNAL = true,
        const bool FLIP_ALL = true,
        const bool USE_DIFFERENT_FLIP_AXES = false
      );

      //! @brief constructor from a set of flip axis, flip internal and and a flip all boolean
      //! @param FLIP_AXES Set of axes to apply the flip around, (randomly chosen at flip time)
      //! @param FLIP_INTERNAL boolean to flip internally
      //! @param FLIP_ALL boolean to decide whether to flip all strands or only few
      //! @param USE_DIFFERENT_FLIP_AXES boolean to decide whether to use different flip axes for each flip
      MutateDomainFlip
      (
        const storage::Set< coord::Axis> &FLIP_AXES,
        const bool FLIP_INTERNAL = true,
        const bool FLIP_ALL = true,
        const bool USE_DIFFERENT_FLIP_AXES = false
      );

      //! @brief Clone function
      //! @return pointer to new MutateDomainFlip
      MutateDomainFlip *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns flip internal or not
      //! @return flip internal or not
      bool GetFlipInternal() const
      {
        return m_FlipInternal;
      }

      //! @brief returns flip all or not
      //! @return flip all or not
      bool GetFlipAll() const
      {
        return m_FlipAll;
      }

      //! @brief return flip axes
      //! @return flip axes
      const storage::Set< coord::Axis> &GetFlipAxes() const;

      //! @brief returns use different flip axes or not
      //! @return use different flip axes or not
      bool GetUseDifferentFlipAxes() const
      {
        return m_UseDifferentFlipAxes;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes a Domain and return a mutated Domain
      //! @param THIS_DOMAIN Domain which will be mutated
      //! @return MutateResult with the mutated Domain
      math::MutateResult< assemble::Domain> operator()( const assemble::Domain &THIS_DOMAIN) const;

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

    }; // class MutateDomainFlip

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_DOMAIN_FLIP_H_ 
