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

#ifndef BCL_ASSEMBLE_BIOMOLECULE_H_
#define BCL_ASSEMBLE_BIOMOLECULE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "quality/bcl_quality_superimpose_interface.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Biomolecule
    //! @brief helps to define a bio molecule
    //! @details Given a protein model, the operator()() finds transformation matrices that will superimpose chains onto
    //!          a single chain by a given superimposer and threshold for a given set of atom coordinates
    //!
    //! @see @link example_assemble_biomolecule.cpp @endlink
    //! @author woetzen
    //! @date Jul 26, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Biomolecule :
      public math::MutateInterface< ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! atom types to use to acquire coordinates
      storage::Set< biol::AtomType> m_AtomTypes;

      //! superimpose measure to use
      util::ShPtr< quality::SuperimposeInterface> m_Superimpose;

      //! threshold as criteria for an identical chain
      double m_Threshold;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Biomolecule();

      //! @brief constructor from required members
      //! @param ATOM_TYPES atoms types to use for superimposition
      //! @param SUPERIMPOSE the superimposition measure to use
      //! @param THRESHOLD the threshold for which two chains are considered the same
      Biomolecule
      (
        const storage::Set< biol::AtomType> &ATOM_TYPES,
        const util::ShPtr< quality::SuperimposeInterface> &SUPERIMPOSE,
        const double THRESHOLD
      );

      //! @brief Clone function
      //! @return pointer to new Biomolecule
      Biomolecule *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a PROTEIN_MODEL and returning a protein model with a protein model multiplier in the
      //! data, that defines a symmetric subunits
      //! @param PROTEIN_MODEL model with multiple chains
      //! @return MutateResult that results from finding the bio molecules
      math::MutateResult< ProteinModel> operator()( const ProteinModel &PROTEIN_MODEL) const;

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

    }; // class Biomolecule

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_BIOMOLECULE_H_
