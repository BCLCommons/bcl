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

#ifndef BCL_FOLD_LOOP_DOMAIN_C_TO_N_H_
#define BCL_FOLD_LOOP_DOMAIN_C_TO_N_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_loop_domain.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopDomainCToN
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_fold_loop_domain_c_to_n.cpp @endlink
    //! @author alexanns
    //! @date Jul 6, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoopDomainCToN :
      public LoopDomain
    {

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoopDomainCToN();

      //! @brief constructor taking parameters to create member variables
      //! @param SEGMENTS list of LoopSegments which will be used to create "m_LoopSegments"
      LoopDomainCToN
      (
        const storage::List< LoopSegment> &SEGMENTS
      );

      //! @brief Clone function
      //! @return pointer to new LoopDomainCToN
      LoopDomainCToN *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the directionality of the loop domain
      //! @return sequence direction in this case e_NTerminal since changes are propagated nterminally
      biol::AASequenceFlexibility::SequenceDirection GetSequenceDirection() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the residues in the loop domain in sequence order including anchor sse anchor residue and the
      //!        created psuedo residue
      //! @return map with the amino acids and a bool indicating whether or not the aa is part of a rigid segment or not
      storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > GetResidues() const;

      //! @brief gives the residue that is attached to the anchor sse
      //!        i.e. the residue that is closest to point of attachment
      //! @return ShPtr to residue that is attached to the anchor sse
      const util::ShPtr< biol::AABase> &GetMostProximalLoopSegmentAA() const;

      //! @brief gives the loop segment that is most distant in sequence to attachment to the anchor sse
      //!        this is the sse that the pseudo residue attaches to
      //! @return sse that the pseudo residue attaches to and is the most distant in sequence from anchor sse
      const assemble::SSE &GetMostDistalLoopSegment() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief sets the phi of a residue in the loop domain
      //! @param AA the amino acid that will be found and mutated
      //! @param PHI the phi angle that residue will be set to
      void SetPhi( const biol::AABase &AA, const double PHI);

      //! @brief sets the psi of a residue in the loop domain
      //! @param AA the amino acid that will be found and mutated
      //! @param PSI the psi angle that residue will be set to
      void SetPsi( const biol::AABase &AA, const double PSI);

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

    }; // class LoopDomainCToN

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOOP_DOMAIN_C_TO_N_H_ 
