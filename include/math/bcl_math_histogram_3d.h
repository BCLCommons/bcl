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

#ifndef BCL_MATH_HISTOGRAM_3D_H_
#define BCL_MATH_HISTOGRAM_3D_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_tensor.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Histogram3D
    //! @brief This is a Histogram3D class to calculate histograms from a storagevector of Pairs of double values
    //!
    //! @see @link example_math_histogram_3d.cpp @endlink
    //! @author mendenjl
    //! @date Jan 05, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Histogram3D :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      storage::VectorND< 2, double> m_LowerUpperBoundariesX; //!< min and max value in the histogram in x direction
      storage::VectorND< 2, double> m_LowerUpperBoundariesY; //!< min and max value in the histogram in y direction
      storage::VectorND< 2, double> m_LowerUpperBoundariesZ; //!< min and max value in the histogram in z direction
      storage::VectorND< 3, double> m_BinSizeXYZ;            //!< GetSize of one bin in the histogram in x, y, z direction
      Tensor< double>               m_Histogram;             //!< calculated histogram
      linal::Vector< double>        m_BinningX;              //!< Binning in x-dimension. Cached for speed
      linal::Vector< double>        m_BinningY;              //!< Binning in y-dimension. Cached for speed
      linal::Vector< double>        m_BinningZ;              //!< Binning in z-dimension. Cached for speed

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
      Histogram3D();

      //! construct Histogram from X and Y and Z starting values MIN_X_Y_Z, from triplet BINSIZE_X_Y_Z and from NUMBER_OF_BINS_X_Y_Z
      //! @brief INITIAL_VAL initial value of all bins; helpful for pseudocount usage
      Histogram3D
      (
        const storage::VectorND< 3, double> &MIN_X_Y_Z,
        const storage::VectorND< 3, double> &BINSIZE_X_Y_Z,
        const storage::VectorND< 3, size_t> &NUMBER_OF_BINS_X_Y_Z,
        const double &INITIAL_VAL = 0.0
      );

      //! copy constructor
      Histogram3D *Clone() const;

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

      //! get number of bins Y
      size_t GetNumberOfBinsZ() const;

      //! get binsize XY
      const storage::VectorND< 3, double> &GetBinSizeXYZ() const;

      //! GetBoundaries X
      const storage::VectorND< 2, double> &GetBoundariesX() const;

      //! GetBoundaries Y
      const storage::VectorND< 2, double> &GetBoundariesY() const;

      //! GetBoundaries Z
      const storage::VectorND< 2, double> &GetBoundariesZ() const;

      //! GetHistogram
      const Tensor< double> &GetHistogram() const;

      //! GetHistogram
      Tensor< double> &GetChangeableHistogram();

      //! get sum of all counts
      double GetSumOfAllCounts() const;

      //! GetBinning
      storage::VectorND< 3, linal::Vector< double> > GetBinningXYZ() const;

    ////////////////
    // operations //
    ////////////////

      //! reset all counts
      void Reset();

      //! pushback a pair of values to the right position in the histogram
      void PushBack( const double &X, const double &Y, const double &Z, const double &WEIGHT = double( 1.0));

      //! Use tri-linear interpolation to obtain a value at an arbitrary point
      //! @param X, Y, Z the coordinates of interest
      double Interpolate( const double &X, const double &Y, const double &Z) const;

      //! Get the value for the nearest bin
      //! @param X, Y, Z the coordinates of interest
      double Value( const double &X, const double &Y, const double &Z) const;

      //! @brief combine this with a given histogram2d by adding up all counts
      //! all parameters have to be identical for this operation to work
      //! @param HISTOGRAM_2D histogram to be added to this one
      //! @return true if it was successful
      bool Combine( const Histogram3D &HISTOGRAM_2D);

      //! checks if there is any count
      bool IsEmpty() const;

      //! @brief normalize the counts in the histogram
      //! each bin is divided by the total number of counts (even the boundary counts)
      void Normalize();

      //! @brief Normalize by another histogram 2d, which represents background counts
      //! @param BACKGROUND background counts histogram
      void NormalizeByBackground( const Histogram3D &HIST);

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

      //! read Histogram3D from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief setup bins according to the boundaries
      void SetupBinning();

    }; // class Histogram3D

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_HISTOGRAM_3D_H_
