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

#ifndef BCL_MATH_HISTOGRAM_2D_H_
#define BCL_MATH_HISTOGRAM_2D_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Histogram2D
    //! @brief This is a Histogram2D class to calculate histograms from a storagevector of Pairs of double values
    //!
    //! @see @link example_math_histogram_2d.cpp @endlink
    //! @author woetzen
    //! @date 11.05.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Histogram2D :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      storage::VectorND< 2, double> m_LowerUpperBoundariesX; //!< min and max value in the histogram in x direction
      storage::VectorND< 2, double> m_LowerUpperBoundariesY; //!< min and max value in the histogram in y direction
      storage::VectorND< 2, double> m_BinSizeXY;             //!< GetSize of one bin in the histogram in x and y direction
      linal::Matrix< double>        m_Histogram;             //!< calculated histogram

      //! @brief string to indicate the center
      //! @return string to indicate the center
      static const std::string &GetCenterString();

      //! @brief string to indicate counts
      //! @return string to indicate counts
      static const std::string &GetCountsString();

      //! @brief string to indicate the left boundary in output
      //! @return string to indicate the left boundary in output
      static const std::string &GetLeftBoundaryString();

      //! @brief string to indicate a bin
      //! @return string to indicate a bin
      static const std::string &GetBinString();

      //! @brief string to indicate the right boundary in output
      //! @return string to indicate the right boundary in output
      static const std::string &GetRightBoundaryString();

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor - initializes everything to 0
      Histogram2D();

      //! construct Histogram from X and Y starting values MIN_X_Y, from Pair BINSIZE_X_Y and from NUMBER_OF_BINS_X_Y
      Histogram2D
      (
        const storage::VectorND< 2, double> &MIN_X_Y,
        const storage::VectorND< 2, double> &BINSIZE_X_Y,
        const storage::VectorND< 2, size_t> &NUMBER_OF_BINS_X_Y
      );

      //! copy constructor
      Histogram2D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! get number of bins X
      size_t GetNumberOfBinsX() const;

      //! get number of bins Y
      size_t GetNumberOfBinsY() const;

      //! get binsize XY
      const storage::VectorND< 2, double> &GetBinSizeXY() const;

      //! GetBoundaries X
      const storage::VectorND< 2, double> &GetBoundariesX() const;

      //! GetBoundaries Y
      const storage::VectorND< 2, double> &GetBoundariesY() const;

      //! GetBoundaries counts X
      storage::VectorND< 2, linal::Vector< double> > GetBoundariesCountsX() const;

      //! GetBoundaries counts Y
      storage::VectorND< 2, linal::Vector< double> > GetBoundariesCountsY() const;

      //! GetHistogram
      linal::Matrix< double> GetHistogram() const;

      //! get sum of all counts
      double GetSumOfAllCounts() const;

      //! GetBinning
      storage::VectorND< 2, linal::Vector< double> > GetBinningXY() const;

    ////////////////
    // operations //
    ////////////////

      //! calculate the Histogram2D from a VALUES_VECTOR
      Histogram2D &CalculateHistogram( const storage::Vector< storage::VectorND< 2, double> > &VALUES_VECTOR);

      //! reset all counts
      void Reset();

      //! pushback a pair of values to the right position in the histogram
      void PushBack( const storage::VectorND< 2, double> &PAIR_OF_VALUES, const double &WEIGHT = double( 1.0));

      //! @brief combine this with a given histogram2d by adding up all counts
      //! all parameters have to be identical for this operation to work
      //! @param HISTOGRAM_2D histogram to be added to this one
      //! @return true if it was successful
      bool Combine( const Histogram2D &HISTOGRAM_2D);

      //! checks if there is any count
      bool IsEmpty() const;

      //! @brief normalize the counts in the histogram
      //! each bin is divided by the total number of counts (even the boundary counts)
      void Normalize();

      //! @brief normalize the row counts in the histogram
      //! @detail each cell is divided by the maximum number of counts in the row
      void NormalizeRows();

      //! @brief normalize the column counts in the histogram
      //! each cell is divided by the total number of counts in the column
      void NormalizeY();

      //! @brief Normalize by another histogram 2d, which represents background counts
      //! @param BACKGROUND background counts histogram
      void NormalizeByBackground( const Histogram2D &HIST);

      //! @brief calculate the logarithm of the column counts + 1 in the histogram
      void Log();

    //////////////////////
    // input and output //
    //////////////////////

      //! write into std::ostream x and y values horizontally
      std::ostream &Write( std::ostream &OSTREAM, const util::Format &FORMAT_BINNING, const util::Format &FORMAT_VALUES) const;

    protected:

      //! write into std::ostream x and y values horizontally
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return Write( OSTREAM, util::Format().W( 8).FFP( 3), util::Format().W( 8).FFP( 3));
      }

      //! read Histogram2D from std::istream
      std::istream &Read( std::istream &ISTREAM);

    }; // class Histogram2D

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_HISTOGRAM_2D_H_
