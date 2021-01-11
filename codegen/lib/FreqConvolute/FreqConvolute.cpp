//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  FreqConvolute.cpp
//
//  Code generation for function 'FreqConvolute'
//


// Include files
#include "FreqConvolute.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "rt_nonfinite.h"
#include <cmath>

// Type Definitions
namespace coder
{
  namespace internal
  {
    class FFTImplementationCallback
    {
     public:
      static void get_algo_sizes(int nfft, boolean_T useRadix2, int *n2blue, int
        *nRows);
      static void dobluesteinfft(const ::coder::array<double, 1U> &x, int n2blue,
        int nfft, const ::coder::array<double, 2U> &costab, const ::coder::array<
        double, 2U> &sintab, const ::coder::array<double, 2U> &sintabinv, ::
        coder::array<creal_T, 1U> &y);
      static void r2br_r2dit_trig_impl(const ::coder::array<creal_T, 1U> &x, int
        unsigned_nRows, const ::coder::array<double, 2U> &costab, const ::coder::
        array<double, 2U> &sintab, ::coder::array<creal_T, 1U> &y);
      static void doHalfLengthRadix2(const ::coder::array<double, 1U> &x, ::
        coder::array<creal_T, 1U> &y, int unsigned_nRows, const ::coder::array<
        double, 2U> &costab, const ::coder::array<double, 2U> &sintab);
     protected:
      static void doHalfLengthBluestein(const ::coder::array<double, 1U> &x, ::
        coder::array<creal_T, 1U> &y, int nrowsx, int nRows, int nfft, const ::
        coder::array<creal_T, 1U> &wwc, const ::coder::array<double, 2U> &costab,
        const ::coder::array<double, 2U> &sintab, const ::coder::array<double,
        2U> &costabinv, const ::coder::array<double, 2U> &sintabinv);
    };
  }
}

// Function Declarations
namespace coder
{
  static void fft(const ::coder::array<double, 1U> &x, double varargin_1, ::
                  coder::array<creal_T, 1U> &y);
}

static double rt_powd_snf(double u0, double u1);

// Function Definitions
namespace coder
{
  namespace internal
  {
    void FFTImplementationCallback::doHalfLengthBluestein(const ::coder::array<
      double, 1U> &x, ::coder::array<creal_T, 1U> &y, int nrowsx, int nRows, int
      nfft, const ::coder::array<creal_T, 1U> &wwc, const ::coder::array<double,
      2U> &costab, const ::coder::array<double, 2U> &sintab, const ::coder::
      array<double, 2U> &costabinv, const ::coder::array<double, 2U> &sintabinv)
    {
      array<creal_T, 1U> fv;
      array<creal_T, 1U> fy;
      array<creal_T, 1U> reconVar1;
      array<creal_T, 1U> reconVar2;
      array<creal_T, 1U> ytmp;
      array<double, 2U> b_costab;
      array<double, 2U> b_sintab;
      array<double, 2U> costab1q;
      array<double, 2U> hcostabinv;
      array<double, 2U> hsintab;
      array<double, 2U> hsintabinv;
      array<int, 2U> wrapIndex;
      double e;
      double temp_im;
      double temp_re;
      double twid_im;
      double twid_re;
      double z;
      int hnRows;
      int i;
      int iDelta2;
      int ihi;
      int ix;
      int j;
      int ju;
      int k;
      int nRowsD2;
      int nd2;
      int temp_re_tmp;
      boolean_T tst;
      hnRows = nRows / 2;
      ytmp.set_size(hnRows);
      if (hnRows > nrowsx) {
        ytmp.set_size(hnRows);
        for (iDelta2 = 0; iDelta2 < hnRows; iDelta2++) {
          ytmp[iDelta2].re = 0.0;
          ytmp[iDelta2].im = 0.0;
        }
      }

      if ((x.size(0) & 1) == 0) {
        tst = true;
        ju = x.size(0);
      } else if (x.size(0) >= nRows) {
        tst = true;
        ju = nRows;
      } else {
        tst = false;
        ju = x.size(0) - 1;
      }

      if (ju >= nRows) {
        ju = nRows;
      }

      nd2 = nRows << 1;
      e = 6.2831853071795862 / static_cast<double>(nd2);
      ihi = nd2 / 2 / 2;
      costab1q.set_size(1, (ihi + 1));
      costab1q[0] = 1.0;
      nd2 = ihi / 2 - 1;
      for (k = 0; k <= nd2; k++) {
        costab1q[k + 1] = std::cos(e * (static_cast<double>(k) + 1.0));
      }

      iDelta2 = nd2 + 2;
      nd2 = ihi - 1;
      for (k = iDelta2; k <= nd2; k++) {
        costab1q[k] = std::sin(e * static_cast<double>(ihi - k));
      }

      costab1q[ihi] = 0.0;
      ihi = costab1q.size(1) - 1;
      nd2 = (costab1q.size(1) - 1) << 1;
      b_costab.set_size(1, (nd2 + 1));
      b_sintab.set_size(1, (nd2 + 1));
      b_costab[0] = 1.0;
      b_sintab[0] = 0.0;
      for (k = 0; k < ihi; k++) {
        b_costab[k + 1] = costab1q[k + 1];
        b_sintab[k + 1] = -costab1q[(ihi - k) - 1];
      }

      iDelta2 = costab1q.size(1);
      for (k = iDelta2; k <= nd2; k++) {
        b_costab[k] = -costab1q[nd2 - k];
        b_sintab[k] = -costab1q[k - ihi];
      }

      nd2 = costab.size(1) / 2;
      costab1q.set_size(1, nd2);
      hsintab.set_size(1, nd2);
      hcostabinv.set_size(1, nd2);
      hsintabinv.set_size(1, nd2);
      for (i = 0; i < nd2; i++) {
        iDelta2 = ((i + 1) << 1) - 2;
        costab1q[i] = costab[iDelta2];
        hsintab[i] = sintab[iDelta2];
        hcostabinv[i] = costabinv[iDelta2];
        hsintabinv[i] = sintabinv[iDelta2];
      }

      reconVar1.set_size(hnRows);
      reconVar2.set_size(hnRows);
      nd2 = 0;
      wrapIndex.set_size(1, hnRows);
      for (i = 0; i < hnRows; i++) {
        reconVar1[i].re = b_sintab[nd2] + 1.0;
        reconVar1[i].im = -b_costab[nd2];
        reconVar2[i].re = 1.0 - b_sintab[nd2];
        reconVar2[i].im = b_costab[nd2];
        nd2 += 2;
        if (i + 1 != 1) {
          wrapIndex[i] = (hnRows - i) + 1;
        } else {
          wrapIndex[0] = 1;
        }
      }

      nd2 = 0;
      e = static_cast<double>(ju) / 2.0;
      iDelta2 = static_cast<int>(e);
      for (ix = 0; ix < iDelta2; ix++) {
        temp_re_tmp = (hnRows + ix) - 1;
        temp_re = wwc[temp_re_tmp].re;
        temp_im = wwc[temp_re_tmp].im;
        twid_re = x[nd2];
        twid_im = x[nd2 + 1];
        ytmp[ix].re = temp_re * twid_re + temp_im * twid_im;
        ytmp[ix].im = temp_re * twid_im - temp_im * twid_re;
        nd2 += 2;
      }

      if (!tst) {
        temp_re_tmp = (hnRows + static_cast<int>(e)) - 1;
        temp_re = wwc[temp_re_tmp].re;
        temp_im = wwc[temp_re_tmp].im;
        twid_re = x[nd2];
        ytmp[static_cast<int>(static_cast<double>(ju) / 2.0)].re = temp_re *
          twid_re + temp_im * 0.0;
        ytmp[static_cast<int>(static_cast<double>(ju) / 2.0)].im = temp_re * 0.0
          - temp_im * twid_re;
        if (static_cast<int>(e) + 2 <= hnRows) {
          iDelta2 = static_cast<int>(static_cast<double>(ju) / 2.0) + 2;
          for (i = iDelta2; i <= hnRows; i++) {
            ytmp[i - 1].re = 0.0;
            ytmp[i - 1].im = 0.0;
          }
        }
      } else {
        if (static_cast<int>(e) + 1 <= hnRows) {
          iDelta2 = static_cast<int>(static_cast<double>(ju) / 2.0) + 1;
          for (i = iDelta2; i <= hnRows; i++) {
            ytmp[i - 1].re = 0.0;
            ytmp[i - 1].im = 0.0;
          }
        }
      }

      z = static_cast<double>(nfft) / 2.0;
      fy.set_size((static_cast<int>(z)));
      if (static_cast<int>(z) > ytmp.size(0)) {
        nd2 = static_cast<int>(z);
        fy.set_size((static_cast<int>(z)));
        for (iDelta2 = 0; iDelta2 < nd2; iDelta2++) {
          fy[iDelta2].re = 0.0;
          fy[iDelta2].im = 0.0;
        }
      }

      ju = ytmp.size(0);
      j = static_cast<int>(z);
      if (ju < j) {
        j = ju;
      }

      iDelta2 = static_cast<int>(z) - 2;
      nRowsD2 = static_cast<int>(z) / 2;
      k = nRowsD2 / 2;
      ix = 0;
      nd2 = 0;
      ju = 0;
      for (i = 0; i <= j - 2; i++) {
        fy[nd2] = ytmp[ix];
        ihi = static_cast<int>(z);
        tst = true;
        while (tst) {
          ihi >>= 1;
          ju ^= ihi;
          tst = ((ju & ihi) == 0);
        }

        nd2 = ju;
        ix++;
      }

      fy[nd2] = ytmp[ix];
      if (static_cast<int>(z) > 1) {
        for (i = 0; i <= iDelta2; i += 2) {
          temp_re = fy[i + 1].re;
          temp_im = fy[i + 1].im;
          e = fy[i].re;
          twid_re = fy[i].im;
          fy[i + 1].re = fy[i].re - fy[i + 1].re;
          fy[i + 1].im = fy[i].im - fy[i + 1].im;
          e += temp_re;
          twid_re += temp_im;
          fy[i].re = e;
          fy[i].im = twid_re;
        }
      }

      nd2 = 2;
      iDelta2 = 4;
      ix = ((k - 1) << 2) + 1;
      while (k > 0) {
        for (i = 0; i < ix; i += iDelta2) {
          temp_re_tmp = i + nd2;
          temp_re = fy[temp_re_tmp].re;
          temp_im = fy[temp_re_tmp].im;
          fy[temp_re_tmp].re = fy[i].re - temp_re;
          fy[temp_re_tmp].im = fy[i].im - temp_im;
          fy[i].re = fy[i].re + temp_re;
          fy[i].im = fy[i].im + temp_im;
        }

        ju = 1;
        for (j = k; j < nRowsD2; j += k) {
          twid_re = costab1q[j];
          twid_im = hsintab[j];
          i = ju;
          ihi = ju + ix;
          while (i < ihi) {
            temp_re_tmp = i + nd2;
            temp_re = twid_re * fy[temp_re_tmp].re - twid_im * fy[temp_re_tmp].
              im;
            temp_im = twid_re * fy[temp_re_tmp].im + twid_im * fy[temp_re_tmp].
              re;
            fy[temp_re_tmp].re = fy[i].re - temp_re;
            fy[temp_re_tmp].im = fy[i].im - temp_im;
            fy[i].re = fy[i].re + temp_re;
            fy[i].im = fy[i].im + temp_im;
            i += iDelta2;
          }

          ju++;
        }

        k /= 2;
        nd2 = iDelta2;
        iDelta2 += iDelta2;
        ix -= nd2;
      }

      FFTImplementationCallback::r2br_r2dit_trig_impl((wwc), (static_cast<int>(z)),
        (costab1q), (hsintab), (fv));
      nd2 = fy.size(0);
      for (iDelta2 = 0; iDelta2 < nd2; iDelta2++) {
        twid_re = fy[iDelta2].re * fv[iDelta2].im + fy[iDelta2].im * fv[iDelta2]
          .re;
        fy[iDelta2].re = fy[iDelta2].re * fv[iDelta2].re - fy[iDelta2].im *
          fv[iDelta2].im;
        fy[iDelta2].im = twid_re;
      }

      FFTImplementationCallback::r2br_r2dit_trig_impl((fy), (static_cast<int>(z)),
        (hcostabinv), (hsintabinv), (fv));
      if (fv.size(0) > 1) {
        e = 1.0 / static_cast<double>(fv.size(0));
        nd2 = fv.size(0);
        for (iDelta2 = 0; iDelta2 < nd2; iDelta2++) {
          fv[iDelta2].re = e * fv[iDelta2].re;
          fv[iDelta2].im = e * fv[iDelta2].im;
        }
      }

      nd2 = 0;
      iDelta2 = wwc.size(0);
      for (k = hnRows; k <= iDelta2; k++) {
        ytmp[nd2].re = wwc[k - 1].re * fv[k - 1].re + wwc[k - 1].im * fv[k - 1].
          im;
        ytmp[nd2].im = wwc[k - 1].re * fv[k - 1].im - wwc[k - 1].im * fv[k - 1].
          re;
        nd2++;
      }

      for (i = 0; i < hnRows; i++) {
        iDelta2 = wrapIndex[i];
        e = ytmp[iDelta2 - 1].re;
        twid_re = -ytmp[iDelta2 - 1].im;
        y[i].re = 0.5 * ((ytmp[i].re * reconVar1[i].re - ytmp[i].im *
                          reconVar1[i].im) + (e * reconVar2[i].re - twid_re *
          reconVar2[i].im));
        y[i].im = 0.5 * ((ytmp[i].re * reconVar1[i].im + ytmp[i].im *
                          reconVar1[i].re) + (e * reconVar2[i].im + twid_re *
          reconVar2[i].re));
        iDelta2 = hnRows + i;
        y[iDelta2].re = 0.5 * ((ytmp[i].re * reconVar2[i].re - ytmp[i].im *
          reconVar2[i].im) + (e * reconVar1[i].re - twid_re * reconVar1[i].im));
        y[iDelta2].im = 0.5 * ((ytmp[i].re * reconVar2[i].im + ytmp[i].im *
          reconVar2[i].re) + (e * reconVar1[i].im + twid_re * reconVar1[i].re));
      }
    }

    void FFTImplementationCallback::doHalfLengthRadix2(const ::coder::array<
      double, 1U> &x, ::coder::array<creal_T, 1U> &y, int unsigned_nRows, const ::
      coder::array<double, 2U> &costab, const ::coder::array<double, 2U> &sintab)
    {
      array<creal_T, 1U> reconVar1;
      array<creal_T, 1U> reconVar2;
      array<double, 2U> hcostab;
      array<double, 2U> hsintab;
      array<int, 2U> wrapIndex;
      array<int, 1U> bitrevIndex;
      double temp2_im;
      double temp2_re;
      double temp_im;
      double temp_re;
      double y_im;
      double y_im_tmp;
      double z_tmp;
      int hszCostab;
      int i;
      int ihi;
      int istart;
      int iy;
      int ju;
      int k;
      int nRows;
      int nRowsD2;
      int nRowsM2;
      boolean_T tst;
      nRows = unsigned_nRows / 2;
      ihi = y.size(0);
      if (ihi >= nRows) {
        ihi = nRows;
      }

      nRowsM2 = nRows - 2;
      nRowsD2 = nRows / 2;
      k = nRowsD2 / 2;
      hszCostab = costab.size(1) / 2;
      hcostab.set_size(1, hszCostab);
      hsintab.set_size(1, hszCostab);
      for (i = 0; i < hszCostab; i++) {
        iy = ((i + 1) << 1) - 2;
        hcostab[i] = costab[iy];
        hsintab[i] = sintab[iy];
      }

      reconVar1.set_size(nRows);
      reconVar2.set_size(nRows);
      wrapIndex.set_size(1, nRows);
      for (i = 0; i < nRows; i++) {
        temp2_im = sintab[i];
        temp2_re = costab[i];
        reconVar1[i].re = temp2_im + 1.0;
        reconVar1[i].im = -temp2_re;
        reconVar2[i].re = 1.0 - temp2_im;
        reconVar2[i].im = temp2_re;
        if (i + 1 != 1) {
          wrapIndex[i] = (nRows - i) + 1;
        } else {
          wrapIndex[0] = 1;
        }
      }

      z_tmp = static_cast<double>(unsigned_nRows) / 2.0;
      ju = 0;
      iy = 1;
      bitrevIndex.set_size((static_cast<int>(z_tmp)));
      hszCostab = static_cast<int>(z_tmp);
      for (istart = 0; istart < hszCostab; istart++) {
        bitrevIndex[istart] = 0;
      }

      for (istart = 0; istart <= ihi - 2; istart++) {
        bitrevIndex[istart] = iy;
        hszCostab = static_cast<int>(z_tmp);
        tst = true;
        while (tst) {
          hszCostab >>= 1;
          ju ^= hszCostab;
          tst = ((ju & hszCostab) == 0);
        }

        iy = ju + 1;
      }

      bitrevIndex[ihi - 1] = iy;
      if ((x.size(0) & 1) == 0) {
        tst = true;
        ihi = x.size(0);
      } else if (x.size(0) >= unsigned_nRows) {
        tst = true;
        ihi = unsigned_nRows;
      } else {
        tst = false;
        ihi = x.size(0) - 1;
      }

      hszCostab = 0;
      if (ihi >= unsigned_nRows) {
        ihi = unsigned_nRows;
      }

      temp2_im = static_cast<double>(ihi) / 2.0;
      istart = static_cast<int>(temp2_im);
      for (i = 0; i < istart; i++) {
        y[bitrevIndex[i] - 1].re = x[hszCostab];
        y[bitrevIndex[i] - 1].im = x[hszCostab + 1];
        hszCostab += 2;
      }

      if (!tst) {
        istart = bitrevIndex[static_cast<int>(temp2_im)] - 1;
        y[istart].re = x[hszCostab];
        y[istart].im = 0.0;
      }

      if (nRows > 1) {
        for (i = 0; i <= nRowsM2; i += 2) {
          temp_re = y[i + 1].re;
          temp_im = y[i + 1].im;
          y[i + 1].re = y[i].re - y[i + 1].re;
          y[i + 1].im = y[i].im - y[i + 1].im;
          y[i].re = y[i].re + temp_re;
          y[i].im = y[i].im + temp_im;
        }
      }

      hszCostab = 2;
      iy = 4;
      ju = ((k - 1) << 2) + 1;
      while (k > 0) {
        for (i = 0; i < ju; i += iy) {
          nRowsM2 = i + hszCostab;
          temp_re = y[nRowsM2].re;
          temp_im = y[nRowsM2].im;
          y[nRowsM2].re = y[i].re - temp_re;
          y[nRowsM2].im = y[i].im - temp_im;
          y[i].re = y[i].re + temp_re;
          y[i].im = y[i].im + temp_im;
        }

        istart = 1;
        for (nRows = k; nRows < nRowsD2; nRows += k) {
          temp2_re = hcostab[nRows];
          temp2_im = hsintab[nRows];
          i = istart;
          ihi = istart + ju;
          while (i < ihi) {
            nRowsM2 = i + hszCostab;
            temp_re = temp2_re * y[nRowsM2].re - temp2_im * y[nRowsM2].im;
            temp_im = temp2_re * y[nRowsM2].im + temp2_im * y[nRowsM2].re;
            y[nRowsM2].re = y[i].re - temp_re;
            y[nRowsM2].im = y[i].im - temp_im;
            y[i].re = y[i].re + temp_re;
            y[i].im = y[i].im + temp_im;
            i += iy;
          }

          istart++;
        }

        k /= 2;
        hszCostab = iy;
        iy += iy;
        ju -= hszCostab;
      }

      hszCostab = static_cast<int>(z_tmp) / 2;
      temp_re = y[0].re;
      temp_im = y[0].im;
      y_im = y[0].re * reconVar1[0].im + y[0].im * reconVar1[0].re;
      temp2_re = y[0].re;
      temp2_im = -y[0].im;
      y[0].re = 0.5 * ((y[0].re * reconVar1[0].re - y[0].im * reconVar1[0].im) +
                       (temp2_re * reconVar2[0].re - temp2_im * reconVar2[0].im));
      y[0].im = 0.5 * (y_im + (temp2_re * reconVar2[0].im + temp2_im *
        reconVar2[0].re));
      y[static_cast<int>(z_tmp)].re = 0.5 * ((temp_re * reconVar2[0].re -
        temp_im * reconVar2[0].im) + (temp_re * reconVar1[0].re - -temp_im *
        reconVar1[0].im));
      y[static_cast<int>(z_tmp)].im = 0.5 * ((temp_re * reconVar2[0].im +
        temp_im * reconVar2[0].re) + (temp_re * reconVar1[0].im + -temp_im *
        reconVar1[0].re));
      for (i = 2; i <= hszCostab; i++) {
        temp_re = y[i - 1].re;
        temp_im = y[i - 1].im;
        istart = wrapIndex[i - 1];
        temp2_re = y[istart - 1].re;
        temp2_im = y[istart - 1].im;
        y_im = y[i - 1].re * reconVar1[i - 1].im + y[i - 1].im * reconVar1[i - 1]
          .re;
        y_im_tmp = -y[istart - 1].im;
        y[i - 1].re = 0.5 * ((y[i - 1].re * reconVar1[i - 1].re - y[i - 1].im *
                              reconVar1[i - 1].im) + (temp2_re * reconVar2[i - 1]
          .re - y_im_tmp * reconVar2[i - 1].im));
        y[i - 1].im = 0.5 * (y_im + (temp2_re * reconVar2[i - 1].im + y_im_tmp *
          reconVar2[i - 1].re));
        iy = (static_cast<int>(z_tmp) + i) - 1;
        y[iy].re = 0.5 * ((temp_re * reconVar2[i - 1].re - temp_im * reconVar2[i
                           - 1].im) + (temp2_re * reconVar1[i - 1].re -
          -temp2_im * reconVar1[i - 1].im));
        y[iy].im = 0.5 * ((temp_re * reconVar2[i - 1].im + temp_im * reconVar2[i
                           - 1].re) + (temp2_re * reconVar1[i - 1].im +
          -temp2_im * reconVar1[i - 1].re));
        y[istart - 1].re = 0.5 * ((temp2_re * reconVar1[istart - 1].re -
          temp2_im * reconVar1[istart - 1].im) + (temp_re * reconVar2[istart - 1]
          .re - -temp_im * reconVar2[istart - 1].im));
        y[istart - 1].im = 0.5 * ((temp2_re * reconVar1[istart - 1].im +
          temp2_im * reconVar1[istart - 1].re) + (temp_re * reconVar2[istart - 1]
          .im + -temp_im * reconVar2[istart - 1].re));
        iy = (istart + static_cast<int>(z_tmp)) - 1;
        y[iy].re = 0.5 * ((temp2_re * reconVar2[istart - 1].re - temp2_im *
                           reconVar2[istart - 1].im) + (temp_re *
          reconVar1[istart - 1].re - -temp_im * reconVar1[istart - 1].im));
        y[iy].im = 0.5 * ((temp2_re * reconVar2[istart - 1].im + temp2_im *
                           reconVar2[istart - 1].re) + (temp_re *
          reconVar1[istart - 1].im + -temp_im * reconVar1[istart - 1].re));
      }

      if (hszCostab != 0) {
        temp2_im = y[hszCostab].re;
        temp_im = y[hszCostab].im;
        y_im = y[hszCostab].re * reconVar1[hszCostab].im + y[hszCostab].im *
          reconVar1[hszCostab].re;
        y_im_tmp = -y[hszCostab].im;
        y[hszCostab].re = 0.5 * ((y[hszCostab].re * reconVar1[hszCostab].re -
          y[hszCostab].im * reconVar1[hszCostab].im) + (temp2_im *
          reconVar2[hszCostab].re - y_im_tmp * reconVar2[hszCostab].im));
        y[hszCostab].im = 0.5 * (y_im + (temp2_im * reconVar2[hszCostab].im +
          y_im_tmp * reconVar2[hszCostab].re));
        istart = static_cast<int>(z_tmp) + hszCostab;
        y[istart].re = 0.5 * ((temp2_im * reconVar2[hszCostab].re - temp_im *
          reconVar2[hszCostab].im) + (temp2_im * reconVar1[hszCostab].re -
          -temp_im * reconVar1[hszCostab].im));
        y[istart].im = 0.5 * ((temp2_im * reconVar2[hszCostab].im + temp_im *
          reconVar2[hszCostab].re) + (temp2_im * reconVar1[hszCostab].im +
          -temp_im * reconVar1[hszCostab].re));
      }
    }

    void FFTImplementationCallback::dobluesteinfft(const ::coder::array<double,
      1U> &x, int n2blue, int nfft, const ::coder::array<double, 2U> &costab,
      const ::coder::array<double, 2U> &sintab, const ::coder::array<double, 2U>
      &sintabinv, ::coder::array<creal_T, 1U> &y)
    {
      array<creal_T, 1U> fv;
      array<creal_T, 1U> fy;
      array<creal_T, 1U> wwc;
      double nt_re;
      int idx;
      int k;
      int nInt2;
      int nInt2m1;
      int rt;
      if ((nfft != 1) && ((nfft & 1) == 0)) {
        int nRows;
        nRows = nfft / 2;
        nInt2m1 = (nRows + nRows) - 1;
        wwc.set_size(nInt2m1);
        idx = nRows;
        rt = 0;
        wwc[nRows - 1].re = 1.0;
        wwc[nRows - 1].im = 0.0;
        nInt2 = nRows << 1;
        for (k = 0; k <= nRows - 2; k++) {
          double nt_im;
          int b_y;
          b_y = ((k + 1) << 1) - 1;
          if (nInt2 - rt <= b_y) {
            rt += b_y - nInt2;
          } else {
            rt += b_y;
          }

          nt_im = -3.1415926535897931 * static_cast<double>(rt) / static_cast<
            double>(nRows);
          if (nt_im == 0.0) {
            nt_re = 1.0;
            nt_im = 0.0;
          } else {
            nt_re = std::cos(nt_im);
            nt_im = std::sin(nt_im);
          }

          wwc[idx - 2].re = nt_re;
          wwc[idx - 2].im = -nt_im;
          idx--;
        }

        idx = 0;
        rt = nInt2m1 - 1;
        for (k = rt; k >= nRows; k--) {
          wwc[k] = wwc[idx];
          idx++;
        }
      } else {
        nInt2m1 = (nfft + nfft) - 1;
        wwc.set_size(nInt2m1);
        idx = nfft;
        rt = 0;
        wwc[nfft - 1].re = 1.0;
        wwc[nfft - 1].im = 0.0;
        nInt2 = nfft << 1;
        for (k = 0; k <= nfft - 2; k++) {
          double nt_im;
          int b_y;
          b_y = ((k + 1) << 1) - 1;
          if (nInt2 - rt <= b_y) {
            rt += b_y - nInt2;
          } else {
            rt += b_y;
          }

          nt_im = -3.1415926535897931 * static_cast<double>(rt) / static_cast<
            double>(nfft);
          if (nt_im == 0.0) {
            nt_re = 1.0;
            nt_im = 0.0;
          } else {
            nt_re = std::cos(nt_im);
            nt_im = std::sin(nt_im);
          }

          wwc[idx - 2].re = nt_re;
          wwc[idx - 2].im = -nt_im;
          idx--;
        }

        idx = 0;
        rt = nInt2m1 - 1;
        for (k = rt; k >= nfft; k--) {
          wwc[k] = wwc[idx];
          idx++;
        }
      }

      y.set_size(nfft);
      if (nfft > x.size(0)) {
        y.set_size(nfft);
        for (rt = 0; rt < nfft; rt++) {
          y[rt].re = 0.0;
          y[rt].im = 0.0;
        }
      }

      if ((n2blue != 1) && ((nfft & 1) == 0)) {
        FFTImplementationCallback::doHalfLengthBluestein((x), (y), (x.size(0)),
          (nfft), (n2blue), (wwc), (costab), (sintab), (costab), (sintabinv));
      } else {
        nInt2m1 = x.size(0);
        if (nfft < nInt2m1) {
          nInt2m1 = nfft;
        }

        rt = 0;
        for (k = 0; k < nInt2m1; k++) {
          nInt2 = (nfft + k) - 1;
          y[k].re = wwc[nInt2].re * x[rt];
          y[k].im = wwc[nInt2].im * -x[rt];
          rt++;
        }

        rt = nInt2m1 + 1;
        for (k = rt; k <= nfft; k++) {
          y[k - 1].re = 0.0;
          y[k - 1].im = 0.0;
        }

        FFTImplementationCallback::r2br_r2dit_trig_impl((y), (n2blue), (costab),
          (sintab), (fy));
        FFTImplementationCallback::r2br_r2dit_trig_impl((wwc), (n2blue), (costab),
          (sintab), (fv));
        nInt2m1 = fy.size(0);
        for (rt = 0; rt < nInt2m1; rt++) {
          nt_re = fy[rt].re * fv[rt].im + fy[rt].im * fv[rt].re;
          fy[rt].re = fy[rt].re * fv[rt].re - fy[rt].im * fv[rt].im;
          fy[rt].im = nt_re;
        }

        FFTImplementationCallback::r2br_r2dit_trig_impl((fy), (n2blue), (costab),
          (sintabinv), (fv));
        if (fv.size(0) > 1) {
          nt_re = 1.0 / static_cast<double>(fv.size(0));
          nInt2m1 = fv.size(0);
          for (rt = 0; rt < nInt2m1; rt++) {
            fv[rt].re = nt_re * fv[rt].re;
            fv[rt].im = nt_re * fv[rt].im;
          }
        }

        idx = 0;
        rt = wwc.size(0);
        for (k = nfft; k <= rt; k++) {
          y[idx].re = wwc[k - 1].re * fv[k - 1].re + wwc[k - 1].im * fv[k - 1].
            im;
          y[idx].im = wwc[k - 1].re * fv[k - 1].im - wwc[k - 1].im * fv[k - 1].
            re;
          idx++;
        }
      }
    }

    void FFTImplementationCallback::get_algo_sizes(int nfft, boolean_T useRadix2,
      int *n2blue, int *nRows)
    {
      *n2blue = 1;
      if (useRadix2) {
        *nRows = nfft;
      } else {
        if (nfft > 0) {
          int n;
          int pmax;
          n = (nfft + nfft) - 1;
          pmax = 31;
          if (n <= 1) {
            pmax = 0;
          } else {
            int pmin;
            boolean_T exitg1;
            pmin = 0;
            exitg1 = false;
            while ((!exitg1) && (pmax - pmin > 1)) {
              int k;
              int pow2p;
              k = (pmin + pmax) >> 1;
              pow2p = 1 << k;
              if (pow2p == n) {
                pmax = k;
                exitg1 = true;
              } else if (pow2p > n) {
                pmax = k;
              } else {
                pmin = k;
              }
            }
          }

          *n2blue = 1 << pmax;
        }

        *nRows = *n2blue;
      }
    }

    void FFTImplementationCallback::r2br_r2dit_trig_impl(const ::coder::array<
      creal_T, 1U> &x, int unsigned_nRows, const ::coder::array<double, 2U>
      &costab, const ::coder::array<double, 2U> &sintab, ::coder::array<creal_T,
      1U> &y)
    {
      double temp_im;
      double temp_re;
      double twid_im;
      double twid_re;
      int i;
      int iDelta2;
      int iheight;
      int ix;
      int iy;
      int ju;
      int k;
      int nRowsD2;
      y.set_size(unsigned_nRows);
      if (unsigned_nRows > x.size(0)) {
        y.set_size(unsigned_nRows);
        for (iy = 0; iy < unsigned_nRows; iy++) {
          y[iy].re = 0.0;
          y[iy].im = 0.0;
        }
      }

      iDelta2 = x.size(0);
      if (iDelta2 >= unsigned_nRows) {
        iDelta2 = unsigned_nRows;
      }

      iheight = unsigned_nRows - 2;
      nRowsD2 = unsigned_nRows / 2;
      k = nRowsD2 / 2;
      ix = 0;
      iy = 0;
      ju = 0;
      for (i = 0; i <= iDelta2 - 2; i++) {
        boolean_T tst;
        y[iy] = x[ix];
        iy = unsigned_nRows;
        tst = true;
        while (tst) {
          iy >>= 1;
          ju ^= iy;
          tst = ((ju & iy) == 0);
        }

        iy = ju;
        ix++;
      }

      y[iy] = x[ix];
      if (unsigned_nRows > 1) {
        for (i = 0; i <= iheight; i += 2) {
          temp_re = y[i + 1].re;
          temp_im = y[i + 1].im;
          twid_re = y[i].re;
          twid_im = y[i].im;
          y[i + 1].re = y[i].re - y[i + 1].re;
          y[i + 1].im = y[i].im - y[i + 1].im;
          twid_re += temp_re;
          twid_im += temp_im;
          y[i].re = twid_re;
          y[i].im = twid_im;
        }
      }

      iy = 2;
      iDelta2 = 4;
      iheight = ((k - 1) << 2) + 1;
      while (k > 0) {
        int temp_re_tmp;
        for (i = 0; i < iheight; i += iDelta2) {
          temp_re_tmp = i + iy;
          temp_re = y[temp_re_tmp].re;
          temp_im = y[temp_re_tmp].im;
          y[temp_re_tmp].re = y[i].re - temp_re;
          y[temp_re_tmp].im = y[i].im - temp_im;
          y[i].re = y[i].re + temp_re;
          y[i].im = y[i].im + temp_im;
        }

        ix = 1;
        for (ju = k; ju < nRowsD2; ju += k) {
          int ihi;
          twid_re = costab[ju];
          twid_im = sintab[ju];
          i = ix;
          ihi = ix + iheight;
          while (i < ihi) {
            temp_re_tmp = i + iy;
            temp_re = twid_re * y[temp_re_tmp].re - twid_im * y[temp_re_tmp].im;
            temp_im = twid_re * y[temp_re_tmp].im + twid_im * y[temp_re_tmp].re;
            y[temp_re_tmp].re = y[i].re - temp_re;
            y[temp_re_tmp].im = y[i].im - temp_im;
            y[i].re = y[i].re + temp_re;
            y[i].im = y[i].im + temp_im;
            i += iDelta2;
          }

          ix++;
        }

        k /= 2;
        iy = iDelta2;
        iDelta2 += iDelta2;
        iheight -= iy;
      }
    }
  }

  static void fft(const ::coder::array<double, 1U> &x, double varargin_1, ::
                  coder::array<creal_T, 1U> &y)
  {
    array<double, 2U> costab;
    array<double, 2U> costab1q;
    array<double, 2U> sintab;
    array<double, 2U> sintabinv;
    double twid_im;
    double twid_re;
    int i;
    int iDelta2;
    int ihi;
    int ix;
    int j;
    int ju;
    int k;
    int nd2;
    if ((x.size(0) == 0) || (0 == static_cast<int>(varargin_1))) {
      y.set_size((static_cast<int>(varargin_1)));
      nd2 = static_cast<int>(varargin_1);
      for (ix = 0; ix < nd2; ix++) {
        y[ix].re = 0.0;
        y[ix].im = 0.0;
      }
    } else {
      boolean_T useRadix2;
      useRadix2 = ((static_cast<int>(varargin_1) > 0) && ((static_cast<int>
        (varargin_1) & (static_cast<int>(varargin_1) - 1)) == 0));
      internal::FFTImplementationCallback::get_algo_sizes((static_cast<int>
        (varargin_1)), (useRadix2), (&iDelta2), (&nd2));
      twid_re = 6.2831853071795862 / static_cast<double>(nd2);
      ihi = nd2 / 2 / 2;
      costab1q.set_size(1, (ihi + 1));
      costab1q[0] = 1.0;
      nd2 = ihi / 2 - 1;
      for (k = 0; k <= nd2; k++) {
        costab1q[k + 1] = std::cos(twid_re * (static_cast<double>(k) + 1.0));
      }

      ix = nd2 + 2;
      nd2 = ihi - 1;
      for (k = ix; k <= nd2; k++) {
        costab1q[k] = std::sin(twid_re * static_cast<double>(ihi - k));
      }

      costab1q[ihi] = 0.0;
      if (!useRadix2) {
        ihi = costab1q.size(1) - 1;
        nd2 = (costab1q.size(1) - 1) << 1;
        costab.set_size(1, (nd2 + 1));
        sintab.set_size(1, (nd2 + 1));
        costab[0] = 1.0;
        sintab[0] = 0.0;
        sintabinv.set_size(1, (nd2 + 1));
        for (k = 0; k < ihi; k++) {
          sintabinv[k + 1] = costab1q[(ihi - k) - 1];
        }

        ix = costab1q.size(1);
        for (k = ix; k <= nd2; k++) {
          sintabinv[k] = costab1q[k - ihi];
        }

        for (k = 0; k < ihi; k++) {
          costab[k + 1] = costab1q[k + 1];
          sintab[k + 1] = -costab1q[(ihi - k) - 1];
        }

        ix = costab1q.size(1);
        for (k = ix; k <= nd2; k++) {
          costab[k] = -costab1q[nd2 - k];
          sintab[k] = -costab1q[k - ihi];
        }
      } else {
        ihi = costab1q.size(1) - 1;
        nd2 = (costab1q.size(1) - 1) << 1;
        costab.set_size(1, (nd2 + 1));
        sintab.set_size(1, (nd2 + 1));
        costab[0] = 1.0;
        sintab[0] = 0.0;
        for (k = 0; k < ihi; k++) {
          costab[k + 1] = costab1q[k + 1];
          sintab[k + 1] = -costab1q[(ihi - k) - 1];
        }

        ix = costab1q.size(1);
        for (k = ix; k <= nd2; k++) {
          costab[k] = -costab1q[nd2 - k];
          sintab[k] = -costab1q[k - ihi];
        }

        sintabinv.set_size(1, 0);
      }

      if (useRadix2) {
        y.set_size((static_cast<int>(varargin_1)));
        if (static_cast<int>(varargin_1) > x.size(0)) {
          nd2 = static_cast<int>(varargin_1);
          y.set_size((static_cast<int>(varargin_1)));
          for (ix = 0; ix < nd2; ix++) {
            y[ix].re = 0.0;
            y[ix].im = 0.0;
          }
        }

        if (static_cast<int>(varargin_1) != 1) {
          internal::FFTImplementationCallback::doHalfLengthRadix2((x), (y), (
            static_cast<int>(varargin_1)), (costab), (sintab));
        } else {
          double temp_im;
          double temp_re;
          int nRowsD2;
          nd2 = x.size(0);
          j = static_cast<int>(varargin_1);
          if (nd2 < j) {
            j = nd2;
          }

          iDelta2 = static_cast<int>(varargin_1) - 2;
          nRowsD2 = static_cast<int>(varargin_1) / 2;
          k = nRowsD2 / 2;
          ix = 0;
          nd2 = 0;
          ju = 0;
          for (i = 0; i <= j - 2; i++) {
            y[nd2].re = x[ix];
            y[nd2].im = 0.0;
            ihi = static_cast<int>(varargin_1);
            useRadix2 = true;
            while (useRadix2) {
              ihi >>= 1;
              ju ^= ihi;
              useRadix2 = ((ju & ihi) == 0);
            }

            nd2 = ju;
            ix++;
          }

          y[nd2].re = x[ix];
          y[nd2].im = 0.0;
          if (static_cast<int>(varargin_1) > 1) {
            for (i = 0; i <= iDelta2; i += 2) {
              temp_re = y[i + 1].re;
              temp_im = y[i + 1].im;
              twid_re = y[i].re;
              twid_im = y[i].im;
              y[i + 1].re = y[i].re - y[i + 1].re;
              y[i + 1].im = y[i].im - y[i + 1].im;
              twid_re += temp_re;
              twid_im += temp_im;
              y[i].re = twid_re;
              y[i].im = twid_im;
            }
          }

          nd2 = 2;
          iDelta2 = 4;
          ix = ((k - 1) << 2) + 1;
          while (k > 0) {
            int temp_re_tmp;
            for (i = 0; i < ix; i += iDelta2) {
              temp_re_tmp = i + nd2;
              temp_re = y[temp_re_tmp].re;
              temp_im = y[temp_re_tmp].im;
              y[temp_re_tmp].re = y[i].re - temp_re;
              y[temp_re_tmp].im = y[i].im - temp_im;
              y[i].re = y[i].re + temp_re;
              y[i].im = y[i].im + temp_im;
            }

            ju = 1;
            for (j = k; j < nRowsD2; j += k) {
              twid_re = costab[j];
              twid_im = sintab[j];
              i = ju;
              ihi = ju + ix;
              while (i < ihi) {
                temp_re_tmp = i + nd2;
                temp_re = twid_re * y[temp_re_tmp].re - twid_im * y[temp_re_tmp]
                  .im;
                temp_im = twid_re * y[temp_re_tmp].im + twid_im * y[temp_re_tmp]
                  .re;
                y[temp_re_tmp].re = y[i].re - temp_re;
                y[temp_re_tmp].im = y[i].im - temp_im;
                y[i].re = y[i].re + temp_re;
                y[i].im = y[i].im + temp_im;
                i += iDelta2;
              }

              ju++;
            }

            k /= 2;
            nd2 = iDelta2;
            iDelta2 += iDelta2;
            ix -= nd2;
          }
        }
      } else {
        internal::FFTImplementationCallback::dobluesteinfft((x), (iDelta2), (
          static_cast<int>(varargin_1)), (costab), (sintab), (sintabinv), (y));
      }
    }
  }
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = std::pow(u0, u1);
    }
  }

  return y;
}

void FreqConvolute(coder::array<double, 1U> &input, const coder::array<double,
                   1U> &ir_frame_real, const coder::array<double, 1U>
                   &ir_frame_imag, double fft_frame_size, coder::array<double,
                   1U> &out, double *out_size)
{
  coder::array<creal_T, 1U> X;
  coder::array<creal_T, 1U> fv;
  coder::array<creal_T, 1U> wwc;
  coder::array<creal_T, 1U> y;
  coder::array<double, 2U> costab;
  coder::array<double, 2U> costab1q;
  coder::array<double, 2U> sintab;
  coder::array<double, 2U> sintabinv;
  coder::array<double, 1U> b_input;
  double im;
  double nt_im;
  double nt_re;
  int N2blue;
  int i;
  int n2;
  int nd2;
  int nfft;
  n2 = static_cast<int>(fft_frame_size - static_cast<double>(input.size(0)));
  b_input.set_size((input.size(0) + n2));
  nd2 = input.size(0);
  for (i = 0; i < nd2; i++) {
    b_input[i] = input[i];
  }

  for (i = 0; i < n2; i++) {
    b_input[i + input.size(0)] = 0.0;
  }

  input.set_size(b_input.size(0));
  n2 = b_input.size(0);
  for (i = 0; i < n2; i++) {
    input[i] = b_input[i];
  }

  coder::fft(input, fft_frame_size, X);
  n2 = X.size(0);
  for (i = 0; i < n2; i++) {
    nt_re = ir_frame_real[i];
    nt_im = ir_frame_imag[i];
    im = X[i].re * nt_im + X[i].im * nt_re;
    X[i].re = X[i].re * nt_re - X[i].im * nt_im;
    X[i].im = im;
  }

  nfft = static_cast<int>(fft_frame_size);
  if ((X.size(0) == 0) || (0 == static_cast<int>(fft_frame_size))) {
    y.set_size((static_cast<int>(fft_frame_size)));
    n2 = static_cast<int>(fft_frame_size);
    for (i = 0; i < n2; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  } else {
    int k;
    int nInt2;
    boolean_T useRadix2;
    useRadix2 = ((static_cast<int>(fft_frame_size) > 0) && ((static_cast<int>
      (fft_frame_size) & (static_cast<int>(fft_frame_size) - 1)) == 0));
    coder::internal::FFTImplementationCallback::get_algo_sizes((static_cast<int>
      (fft_frame_size)), (useRadix2), (&N2blue), (&nd2));
    nt_re = 6.2831853071795862 / static_cast<double>(nd2);
    nInt2 = nd2 / 2 / 2;
    costab1q.set_size(1, (nInt2 + 1));
    costab1q[0] = 1.0;
    nd2 = nInt2 / 2 - 1;
    for (k = 0; k <= nd2; k++) {
      costab1q[k + 1] = std::cos(nt_re * (static_cast<double>(k) + 1.0));
    }

    i = nd2 + 2;
    nd2 = nInt2 - 1;
    for (k = i; k <= nd2; k++) {
      costab1q[k] = std::sin(nt_re * static_cast<double>(nInt2 - k));
    }

    costab1q[nInt2] = 0.0;
    if (!useRadix2) {
      nInt2 = costab1q.size(1) - 1;
      n2 = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, (n2 + 1));
      sintab.set_size(1, (n2 + 1));
      costab[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, (n2 + 1));
      for (k = 0; k < nInt2; k++) {
        sintabinv[k + 1] = costab1q[(nInt2 - k) - 1];
      }

      i = costab1q.size(1);
      for (k = i; k <= n2; k++) {
        sintabinv[k] = costab1q[k - nInt2];
      }

      for (k = 0; k < nInt2; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(nInt2 - k) - 1];
      }

      i = costab1q.size(1);
      for (k = i; k <= n2; k++) {
        costab[k] = -costab1q[n2 - k];
        sintab[k] = -costab1q[k - nInt2];
      }
    } else {
      nInt2 = costab1q.size(1) - 1;
      n2 = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, (n2 + 1));
      sintab.set_size(1, (n2 + 1));
      costab[0] = 1.0;
      sintab[0] = 0.0;
      for (k = 0; k < nInt2; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = costab1q[(nInt2 - k) - 1];
      }

      i = costab1q.size(1);
      for (k = i; k <= n2; k++) {
        costab[k] = -costab1q[n2 - k];
        sintab[k] = costab1q[k - nInt2];
      }

      sintabinv.set_size(1, 0);
    }

    if (useRadix2) {
      coder::internal::FFTImplementationCallback::r2br_r2dit_trig_impl((X), (
        static_cast<int>(fft_frame_size)), (costab), (sintab), (y));
      if (y.size(0) > 1) {
        nt_re = 1.0 / static_cast<double>(y.size(0));
        n2 = y.size(0);
        for (i = 0; i < n2; i++) {
          y[i].re = nt_re * y[i].re;
          y[i].im = nt_re * y[i].im;
        }
      }
    } else {
      int idx;
      nd2 = (static_cast<int>(fft_frame_size) + static_cast<int>(fft_frame_size))
        - 1;
      wwc.set_size(nd2);
      idx = static_cast<int>(fft_frame_size);
      n2 = 0;
      wwc[static_cast<int>(fft_frame_size) - 1].re = 1.0;
      wwc[static_cast<int>(fft_frame_size) - 1].im = 0.0;
      nInt2 = static_cast<int>(fft_frame_size) << 1;
      i = static_cast<int>(fft_frame_size);
      for (k = 0; k <= i - 2; k++) {
        int b_y;
        b_y = ((k + 1) << 1) - 1;
        if (nInt2 - n2 <= b_y) {
          n2 += b_y - nInt2;
        } else {
          n2 += b_y;
        }

        nt_im = 3.1415926535897931 * static_cast<double>(n2) / static_cast<
          double>(static_cast<int>(fft_frame_size));
        if (nt_im == 0.0) {
          nt_re = 1.0;
          nt_im = 0.0;
        } else {
          nt_re = std::cos(nt_im);
          nt_im = std::sin(nt_im);
        }

        wwc[idx - 2].re = nt_re;
        wwc[idx - 2].im = -nt_im;
        idx--;
      }

      idx = 0;
      i = nd2 - 1;
      for (k = i; k >= nfft; k--) {
        wwc[k] = wwc[idx];
        idx++;
      }

      y.set_size((static_cast<int>(fft_frame_size)));
      if (static_cast<int>(fft_frame_size) > X.size(0)) {
        n2 = static_cast<int>(fft_frame_size);
        y.set_size((static_cast<int>(fft_frame_size)));
        for (i = 0; i < n2; i++) {
          y[i].re = 0.0;
          y[i].im = 0.0;
        }
      }

      nd2 = static_cast<int>(fft_frame_size);
      nInt2 = X.size(0);
      if (nd2 < nInt2) {
        nInt2 = nd2;
      }

      nd2 = 0;
      for (k = 0; k < nInt2; k++) {
        n2 = (static_cast<int>(fft_frame_size) + k) - 1;
        nt_re = wwc[n2].re;
        nt_im = wwc[n2].im;
        y[k].re = nt_re * X[nd2].re + nt_im * X[nd2].im;
        y[k].im = nt_re * X[nd2].im - nt_im * X[nd2].re;
        nd2++;
      }

      i = nInt2 + 1;
      for (k = i; k <= nfft; k++) {
        y[k - 1].re = 0.0;
        y[k - 1].im = 0.0;
      }

      coder::internal::FFTImplementationCallback::r2br_r2dit_trig_impl((y),
        (N2blue), (costab), (sintab), (X));
      coder::internal::FFTImplementationCallback::r2br_r2dit_trig_impl((wwc),
        (N2blue), (costab), (sintab), (fv));
      n2 = X.size(0);
      for (i = 0; i < n2; i++) {
        im = X[i].re * fv[i].im + X[i].im * fv[i].re;
        X[i].re = X[i].re * fv[i].re - X[i].im * fv[i].im;
        X[i].im = im;
      }

      coder::internal::FFTImplementationCallback::r2br_r2dit_trig_impl((X),
        (N2blue), (costab), (sintabinv), (fv));
      if (fv.size(0) > 1) {
        nt_re = 1.0 / static_cast<double>(fv.size(0));
        n2 = fv.size(0);
        for (i = 0; i < n2; i++) {
          fv[i].re = nt_re * fv[i].re;
          fv[i].im = nt_re * fv[i].im;
        }
      }

      idx = 0;
      i = static_cast<int>(fft_frame_size);
      nd2 = wwc.size(0);
      for (k = i; k <= nd2; k++) {
        y[idx].re = wwc[k - 1].re * fv[k - 1].re + wwc[k - 1].im * fv[k - 1].im;
        y[idx].im = wwc[k - 1].re * fv[k - 1].im - wwc[k - 1].im * fv[k - 1].re;
        nt_re = y[idx].re;
        nt_im = y[idx].im;
        if (nt_im == 0.0) {
          nt_re /= static_cast<double>(static_cast<int>(fft_frame_size));
          im = 0.0;
        } else if (nt_re == 0.0) {
          nt_re = 0.0;
          im = nt_im / static_cast<double>(static_cast<int>(fft_frame_size));
        } else {
          nt_re /= static_cast<double>(static_cast<int>(fft_frame_size));
          im = nt_im / static_cast<double>(static_cast<int>(fft_frame_size));
        }

        y[idx].re = nt_re;
        y[idx].im = im;
        idx++;
      }
    }
  }

  out.set_size(y.size(0));
  n2 = y.size(0);
  for (i = 0; i < n2; i++) {
    out[i] = y[i].re;
  }

  nd2 = out.size(0);
  *out_size = nd2;
}

void FreqConvolute_initialize()
{
}

void FreqConvolute_terminate()
{
  // (no terminate code required)
}

void GetUnisonPartitionedIRFrames(const coder::array<double, 1U> &IR, double
  frame_size, double buffer_size, coder::array<double, 1U> &ir_frames_real,
  coder::array<double, 1U> &ir_frames_imag, double *n_frames, double
  *output_size)
{
  coder::array<creal_T, 1U> IR_fft_frames;
  coder::array<creal_T, 1U> IR_frame;
  coder::array<double, 2U> IR_frames;
  coder::array<double, 1U> IR_padded;
  double padding_needed;
  int b_loop_ub;
  int c_loop_ub;
  int d_loop_ub;
  int i;
  int i1;
  int i2;
  unsigned int idx_start;
  int loop_ub;
  int padding_needed_idx_0;
  *n_frames = std::ceil(static_cast<double>(IR.size(0)) / buffer_size);
  IR_frames.set_size((static_cast<int>(buffer_size)), (static_cast<int>
    (*n_frames)));
  loop_ub = static_cast<int>(buffer_size) * static_cast<int>(*n_frames);
  for (i = 0; i < loop_ub; i++) {
    IR_frames[i] = 0.0;
  }

  *output_size = frame_size * *n_frames;
  loop_ub = static_cast<int>(*output_size);
  IR_fft_frames.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    IR_fft_frames[i].re = 0.0;
    IR_fft_frames[i].im = 0.0;
  }

  padding_needed = *n_frames * buffer_size - static_cast<double>(IR.size(0));
  if (1 > IR.size(0)) {
    loop_ub = 0;
  } else {
    loop_ub = IR.size(0);
  }

  padding_needed_idx_0 = static_cast<int>(*n_frames * buffer_size - static_cast<
    double>(IR.size(0)));
  IR_padded.set_size((loop_ub + static_cast<int>(padding_needed)));
  for (i = 0; i < loop_ub; i++) {
    IR_padded[i] = IR[i];
  }

  b_loop_ub = static_cast<int>(padding_needed);
  for (i = 0; i < b_loop_ub; i++) {
    IR_padded[i + loop_ub] = 0.0;
  }

  idx_start = 1U;
  padding_needed = buffer_size;
  i = static_cast<int>(*n_frames);
  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    if (idx_start > padding_needed) {
      i1 = 0;
      i2 = 0;
    } else {
      i1 = static_cast<int>(idx_start) - 1;
      i2 = static_cast<int>(padding_needed);
    }

    loop_ub = i2 - i1;
    for (i2 = 0; i2 < loop_ub; i2++) {
      IR_frames[i2 + IR_frames.size(0) * b_loop_ub] = IR_padded[i1 + i2];
    }

    idx_start = static_cast<unsigned int>(padding_needed) + 1U;
    padding_needed += buffer_size;
  }

  i = static_cast<int>(*n_frames);
  if (0 <= static_cast<int>(*n_frames) - 1) {
    if (1 > IR_frames.size(0)) {
      c_loop_ub = 0;
    } else {
      c_loop_ub = IR_frames.size(0);
    }

    if (1 > IR_frames.size(0)) {
      i1 = 0;
    } else {
      i1 = IR_frames.size(0);
    }

    padding_needed_idx_0 = static_cast<int>(frame_size - static_cast<double>(i1));
    d_loop_ub = padding_needed_idx_0;
  }

  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    IR_padded.set_size((c_loop_ub + padding_needed_idx_0));
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      IR_padded[i1] = IR_frames[i1 + IR_frames.size(0) * b_loop_ub];
    }

    for (i1 = 0; i1 < d_loop_ub; i1++) {
      IR_padded[i1 + c_loop_ub] = 0.0;
    }

    coder::fft(IR_padded, frame_size, IR_frame);
    padding_needed = ((static_cast<double>(b_loop_ub) + 1.0) - 1.0) * frame_size
      + 1.0;
    if (padding_needed > (static_cast<double>(b_loop_ub) + 1.0) * frame_size) {
      i1 = 1;
    } else {
      i1 = static_cast<int>(padding_needed);
    }

    loop_ub = IR_frame.size(0);
    for (i2 = 0; i2 < loop_ub; i2++) {
      IR_fft_frames[(i1 + i2) - 1] = IR_frame[i2];
    }
  }

  ir_frames_real.set_size(IR_fft_frames.size(0));
  loop_ub = IR_fft_frames.size(0);
  for (i = 0; i < loop_ub; i++) {
    ir_frames_real[i] = IR_fft_frames[i].re;
  }

  ir_frames_imag.set_size(IR_fft_frames.size(0));
  loop_ub = IR_fft_frames.size(0);
  for (i = 0; i < loop_ub; i++) {
    ir_frames_imag[i] = IR_fft_frames[i].im;
  }
}

void RemoveTailBelowThreshold(coder::array<double, 1U> &IR, double threshold_dB)
{
  double threshold_mag;
  int b_i;
  int i;
  boolean_T exitg1;

  //  Removes trailing values below threshold
  threshold_mag = rt_powd_snf(10.0, threshold_dB / 20.0);
  i = static_cast<int>(((-1.0 - static_cast<double>(IR.size(0))) + 1.0) / -1.0);
  b_i = 0;
  exitg1 = false;
  while ((!exitg1) && (b_i <= i - 1)) {
    double c_i;
    c_i = static_cast<double>(IR.size(0)) + -static_cast<double>(b_i);
    if (std::abs(IR[static_cast<int>(c_i) - 1]) > threshold_mag) {
      IR.set_size((static_cast<int>(c_i)));
      exitg1 = true;
    } else {
      b_i++;
    }
  }
}

// End of code generation (FreqConvolute.cpp)
