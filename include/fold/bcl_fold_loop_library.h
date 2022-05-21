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

#ifndef BCL_FOLD_LOOP_LIBRARY_H_
#define BCL_FOLD_LOOP_LIBRARY_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_loop_parameters.h"
#include "storage/bcl_storage_hash_map.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically
#include <iomanip>
#include <sstream>

namespace bcl
{
  namespace fold
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopLibrary
    //! @brief Reads in a given library of known loop conformations. The essential parameters of the conformations
    //! are discretized and used as hash key in the created hash table.
    //!
    //! @see @link example_fold_loop_library.cpp @endlink
    //! @author fischea
    //! @date Dec 18, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LoopLibrary :
      public util::SerializableInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! library of loop templates in a hash map accessible by sequence length
      storage::HashMap< size_t, util::ShPtrVector< LoopParameters> > m_LoopTemplatesLength;

    private:

      //! width of the bins to discretize the euclidean distance
      double m_DistanceBinWidth;

      //! width of the bins to discretize the euler angles
      double m_AngleBinWidth;

      //! library of loop templates in a hash map for quick access to matching templates
      storage::HashMap< std::string, util::ShPtrVector< LoopParameters> > m_LoopTemplates;

      //! path to the loop template library file
      std::string m_LibraryFileName;

      //! map of loop libraries
      static storage::HashMap< std::string, util::ShPtr< LoopLibrary> > s_LoopLibraries;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoopLibrary();

      //! @brief construct from loop template library
      //! @param LIBRARY_FILE_NAME path to the loop template library
      //! @param DISTANCE_BIN_WIDTH width of the bins to discretize the euclidean distance
      //! @param ANGLE_BIN_WIDTH width of the bins to discretize the Euler angles
      LoopLibrary
      (
        const std::string &LIBRARY_FILE_NAME,
        double DISTANCE_BIN_WIDTH = GetDefaultDistanceBinWidth(),
        double ANGLE_BIN_WIDTH = GetDefaultAngleBinWidth()
      );

    public:

      //! @brief returns a loop library from the given file
      //! @param LIBRARY_FILE_NAME path to the loop template library
      //! @param DISTANCE_BIN_WIDTH width of the bins to discretize the euclidean distance
      //! @param ANGLE_BIN_WIDTH width of the bins to discretize the Euler angles
      static util::ShPtr< LoopLibrary> CreateLoopLibrary
      (
        const std::string &LIBRARY_FILE_NAME,
        double DISTANCE_BIN_WIDTH = GetDefaultDistanceBinWidth(),
        double ANGLE_BIN_WIDTH = GetDefaultAngleBinWidth()
      );

      //! @brief clone function
      //! @return pointer to a new LoopLibrary
      LoopLibrary *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the default distance bin width
      //! @return default distance bin width
      static double GetDefaultDistanceBinWidth();

      //! @brief returns the default angle bin width
      //! @return default angle bin width
      static double GetDefaultAngleBinWidth();

      //! @brief returns the loop templates in the library
      //! @return the loop templates in the library
      const storage::HashMap< std::string, util::ShPtrVector< LoopParameters> > &GetTemplates() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief read a loop template library from the given input stream
      //! @param ISTREAM input stream from which to read the loop template library
      void ReadLibrary( std::istream &ISTREAM);

      //! @brief returns suitable templates for the given loop parameters
      //! @param LOOP_PARAMS loop parameters for which to find templates
      //! @return loop templates for the given parameters
      util::ShPtrVector< LoopParameters> FindTemplates( const LoopParameters &LOOP_PARAMS) const;

      //! @brief returns suitable templates for the given loop length
      //! @param LENGTH sequence length of the loop
      //! @return loop templates for the given sequence length
      util::ShPtrVector< LoopParameters> FindTemplates( size_t LENGTH) const;

      //! @brief extends the template library by recombining template to form larger templates
      //! @param map in the format <sequence length>:<number of new templates> to specify how many additional
      //! templates should be generated
      void RecombineTemplates( const storage::HashMap< size_t, size_t> &COUNTS);

      //! @brief combines the two given templates to one larger template
      //! @param TEMPLATE_N the n-terminal template to be combined
      //! @param TEMPLATE_C the c-terminal template to be combined
      //! @return the template resulting from combination of the two given templates
      static util::ShPtr< LoopParameters> Combine
      (
        const LoopParameters &TEMPLATE_N, const LoopParameters &TEMPLATE_C
      );

      //! @brief creates a given number of loop templates of the given length
      //! @param LENGTH sequence length of the generated loop templates
      //! @param COUNT how many loop templates to generate
      //! @return generated loop templates
      util::ShPtrVector< LoopParameters> GenerateTemplates( size_t LENGTH, size_t COUNT) const;

      //! @brief estimates how many templates of different lengths are needed to close the loops for the given
      //! protein models
      //! @param MODELS models for which the loop regions shall be constructed
      //! @return estimation of the number of additional templates required
      static util::ShPtr< storage::HashMap< size_t, size_t> > EstimateTemplateRequirement
      (
        const util::ShPtrVector< assemble::ProteinModel> &MODELS
      );

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief computes the hash key for the given loop parameters
      //! @param LOOP_PARAMS loop parameters for which to compute the hash key
      //! @return hash key for the given loop parameters
      std::string ComputeKey( const LoopParameters &LOOP_PARAMS) const;

      //! @brief convert the argument into a hash-able string
      //! @param ARG argument to be converted into a hash-able string
      //! @return hash-able string of the argument
      template< typename t_DataType>
      static std::string ToString( t_DataType ARG)
      {
        std::ostringstream convert;
        convert << std::setw( 3) << std::setfill( '0') << ARG;

        return convert.str();
      }

      //! @brief adds the given loop template to the template library
      //! @param LOOP_PARAMS loop which shall be added to the library
      void AddTemplate( const LoopParameters &LOOP_PARAMS);

    }; // class LoopLibrary

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOOP_LIBRARY_H_
