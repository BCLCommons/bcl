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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_COMPRESS_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_COMPRESS_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelCompress
    //! @brief Mutate that compress/decompress a protein model
    //! @brief This class moves all the SSEs in the given ProteinModel to and away from the center of the protein model
    //! How much each SSE is moved towards or away from the center of the model is proportional to its starting distance
    //! to the center and is determined by the given member compression ranges (min and max)
    //!
    //! @see @link example_fold_mutate_protein_model_compress.cpp @endlink
    //! @author karakam
    //! @date 10/05/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelCompress :
      public math::MutateInterface< assemble::ProteinModel>
    {
    private:

    //////////
    // data //
    //////////

      //! Compression factor pair composed of an lower and upper values ( such as 0.97 and 1.01)
      math::Range< double> m_CompressionRange;

      //! scheme
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelCompress();

      //! @brief constructor from a specified compression factor and a scheme
      //! @param COMPRESSION_FACTOR compression factor
      //! @param SCHEME Scheme to be used
      MutateProteinModelCompress
      (
        const double COMPRESSION_FACTOR,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelCompress>()
      );

      //! @brief constructor from a compression factor range
      //! @param COMPRESSION_FACTOR_RANGE compression factor range composed of an lower and an upper value, and a scheme
      //! @param SCHEME Scheme to be used
      MutateProteinModelCompress
      (
        const math::Range< double> &COMPRESSION_FACTOR_RANGE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelCompress>()
      );

      //! @brief clone
      MutateProteinModelCompress *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns compression range
      //! @return compression range
      const math::Range< double> GetCompressionRange() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL ProteinModel which will be mutated
      //! @return MutateResult ProteinModel after mutation
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
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateProteinModelCompress

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_MUTATE_PROTEIN_MODEL_COMPRESS_H_
