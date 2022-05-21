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

#ifndef BCL_NMR_ROSETTA_NOE_HANDLER_H_
#define BCL_NMR_ROSETTA_NOE_HANDLER_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_handler_base.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RosettaNOEHandler
    //! @brief Handler for reading/writing Rosetta-style NOE files
    //! @details Reads in NMR-Rosetta formatted files used in Abinitio Folding which contain NMR restraint data.
    //!          Can also write out these files NOE files can be found here:
    //!          http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/file_constraints.html
    //!
    //! @see @link example_nmr_rosetta_noe_handler.cpp @endlink
    //! @author weinerbe
    //! @date Aug 9, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RosettaNOEHandler :
      public restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> >
    {

    private:

    //////////
    // data //
    //////////

      //! defines the amount to be added to the NOE restraint for defining the upper limit
      static const double s_UpperLimit;

      //! the lower limit for NOEs
      static const double s_LowerLimit;

      //! determines the closest distance within the sequence two residues in a restraint can be for the restraint to
      //! be stored
      size_t m_SequenceDistance;

      //! determines if the restraint's seq ids need to be recalculated based on their relation to the pdb ids
      size_t m_SequenceOffset;

      //! prefix for KB potential files
      std::string m_Prefix;

      //! weight to give to score
      double m_Weight;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from the sequence distance and offset
      //! @param SEQUENCE_DISTANCE size_t which is the smallest distance in sequence two residues can be if a restraint
      //!        is going to be stored
      //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
      //! @param PREFIX prefix (path) for KB potential files
      //! @param SCORE_WEIGHT score weight to use
      RosettaNOEHandler
      (
        const size_t SEQUENCE_DISTANCE = 0,
        const size_t SEQUENCE_OFFSET = 0,
        const std::string PREFIX = "",
        const double &SCORE_WEIGHT = 1.0
      );

      //! @brief Clone function
      //! @return pointer to new RosettaNOEHandler
      RosettaNOEHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns short name for the class used over the command line
      //! @return short name for the class used over the command line
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reads restraints from a stream
      //! @param ISTREAM the stream the restraints will be read from
      //! @return the read in restraints
      util::ShPtrVector< restraint::AtomDistance> ReadRestraints( std::istream &ISTREAM) const;

      //! @brief writes restraints to a stream
      //! @param OSTREAM the stream the restraints will be written to
      //! @param RESTRAINT the restraint that will be written to the stream
      //! @return stream the restraints were written to
      std::ostream &WriteRestraints
      (
        std::ostream &OSTREAM,
        const util::ShPtrVector< restraint::AtomDistance> &RESTRAINT
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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief gets the base atom type (i.e. CB)
      //! @param ATOM_TYPE atom type to be converted
      //! @return base atom type
      static biol::AtomType GetBaseAtomType( const biol::AtomType &ATOM_TYPE);

    }; // class RosettaNOEHandler

  } // namespace nmr
} // namespace bcl

#endif // BCL_NMR_ROSETTA_NOE_HANDLER_H_
