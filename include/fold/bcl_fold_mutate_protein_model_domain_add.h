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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_DOMAIN_ADD_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_DOMAIN_ADD_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "coord/bcl_coord.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelDomainAdd
    //! @brief Adds a domain to an existing protein model
    //! @details This mutate class adds a domain to a protein model, using a random rotation and translation.  The
    //!          relative positions of the SSEs within the domain are preserved
    //!
    //! @see @link example_fold_mutate_protein_model_domain_add.cpp @endlink
    //! @author weinerbe
    //! @date Sep 15, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelDomainAdd :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! DomainInterface derived class that has the SSEs to be added to the protein model
      util::ShPtr
      <
        find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::DomainInterface>
      > m_CollectorDomain;

      //! bool whether to fit the domain to a fold template prior to adding
      bool m_FitToFoldTemplate;

      //! scheme
      std::string m_Scheme;

    public:

      //! deviation used for the placement
      static const math::Range< double> s_Deviation;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelDomainAdd();

      //! @brief construct from a domain collector and a scheme
      //! @param DOMAIN_COLLECTOR domain collector to be used to get a domain
      //! @param FIT_TO_FOLD_TEMPLATE bool whether to fit the domain to a fold template prior to adding
      //! @param SCHEME Scheme to be used
      MutateProteinModelDomainAdd
      (
        const find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::DomainInterface> &DOMAIN_COLLECTOR,
        const bool FIT_TO_FOLD_TEMPLATE = false,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelDomainAdd>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelDomainAdd
      MutateProteinModelDomainAdd *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes a ProteinModel and return a mutated ProteinModel
      //! @param PROTEIN_MODEL ProteinModel which will be mutated
      //! @return MutateResult with the mutated ProteinModel
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

      //! @brief takes a set of coordinates and returns the farthest distance from the beginning of the line segment
      //! @param LINE_SEGMENT line segment that contains the point and direction to be used
      //! @param COORDS coordinates to be checked
      //! @return the farthest distance from the beginning of the line segment
      static double GetMaxDistanceFromCenter
      (
        const coord::LineSegment3D &LINE_SEGMENT,
        const util::SiPtrVector< const linal::Vector3D> &COORDS
      );

    }; // class MutateProteinModelDomainAdd

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_DOMAIN_ADD_H_ 
