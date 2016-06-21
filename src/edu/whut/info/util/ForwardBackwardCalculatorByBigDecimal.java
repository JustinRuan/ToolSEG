package edu.whut.info.util;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

/**
 * Created by Justin on 2015/11/23.
 */
public class ForwardBackwardCalculatorByBigDecimal extends ForwardBackwardCalculator {
    protected final MathContext m_Precision = new MathContext(512, RoundingMode.HALF_EVEN);
    protected BigDecimal probability;
    protected BigDecimal[][] alpha = null;
    protected BigDecimal[][] beta = null;

    public <O extends Observation>
    ForwardBackwardCalculatorByBigDecimal(List<? extends O> oseq, Hmm<O> hmm) {
        this(oseq, hmm, EnumSet.of(Computation.ALPHA));
    }

    public <O extends Observation>
    ForwardBackwardCalculatorByBigDecimal(List<? extends O> oseq,
                                          Hmm<O> hmm, EnumSet<Computation> flags) {

        if (oseq.isEmpty())
            throw new IllegalArgumentException("Invalid empty sequence");

        if (flags.contains(Computation.ALPHA))
            computeAlpha(hmm, oseq);

        if (flags.contains(Computation.BETA))
            computeBeta(hmm, oseq);

        computeProbability(oseq, hmm, flags);
    }

    /* Computes the content of the alpha array */
    protected <O extends Observation> void
    computeAlpha(Hmm<? super O> hmm, List<O> oseq) {
        alpha = new BigDecimal[oseq.size()][hmm.nbStates()];

        for (int i = 0; i < hmm.nbStates(); i++)
            computeAlphaInit(hmm, oseq.get(0), i);

        Iterator<O> seqIterator = oseq.iterator();
        if (seqIterator.hasNext())
            seqIterator.next();

        for (int t = 1; t < oseq.size(); t++) {
            O observation = seqIterator.next();

            for (int i = 0; i < hmm.nbStates(); i++)
                computeAlphaStep(hmm, observation, t, i);
        }
    }


    /* Computes alpha[0][i] */
    protected <O extends Observation> void
    computeAlphaInit(Hmm<? super O> hmm, O o, int i) {
        alpha[0][i] = BigDecimal.valueOf(hmm.getPi(i) * hmm.getOpdf(i).probability(o));
    }


    /* Computes alpha[t][j] (t > 0) */
    protected <O extends Observation> void
    computeAlphaStep(Hmm<? super O> hmm, O o, int t, int j) {
        BigDecimal sum = new BigDecimal(0.);

        for (int i = 0; i < hmm.nbStates(); i++) {
            //sum += alpha[t-1][i] * hmm.getAij(i, j);
            sum = sum.add(alpha[t - 1][i].multiply(BigDecimal.valueOf(hmm.getAij(i, j)), m_Precision), m_Precision);
        }


        //alpha[t][j] = sum * hmm.getOpdf(j).probability(o);
        alpha[t][j] = sum.multiply(BigDecimal.valueOf(hmm.getOpdf(j).probability(o)), m_Precision);
    }

    protected <O extends Observation> void computeBeta(Hmm<? super O> hmm, List<O> oseq) {
        beta = new BigDecimal[oseq.size()][hmm.nbStates()];

        int t;
        for (t = 0; t < hmm.nbStates(); ++t) {
            beta[oseq.size() - 1][t] = BigDecimal.valueOf(1);
        }

        for (t = oseq.size() - 2; t >= 0; --t) {
            for (int i = 0; i < hmm.nbStates(); ++i) {
                computeBetaStep(hmm, oseq.get(t + 1), t, i);
            }
        }

    }

    protected <O extends Observation> void computeBetaStep(Hmm<? super O> hmm, O o, int t, int i) {
        BigDecimal sum = new BigDecimal(0.);
        for (int j = 0; j < hmm.nbStates(); ++j) {
            // sum += beta[t + 1][j] * hmm.getAij(i, j) * hmm.getOpdf(j).probability(o);
            sum = sum.add(beta[t + 1][j].multiply(BigDecimal.valueOf(hmm.getAij(i, j)), m_Precision).multiply(BigDecimal.valueOf(hmm.getOpdf(j).probability(o)), m_Precision), m_Precision);
        }

        this.beta[t][i] = sum;
    }

    private <O extends Observation> void
    computeProbability(List<O> oseq, Hmm<? super O> hmm,
                       EnumSet<Computation> flags) {
        probability = new BigDecimal(0);

        if (flags.contains(Computation.ALPHA))
            for (int i = 0; i < hmm.nbStates(); i++) {
                //probability += alpha[oseq.size()-1][i];
                probability = probability.add(alpha[oseq.size() - 1][i], m_Precision);
            }
//        else
//            for (int i = 0; i < hmm.nbStates(); i++)
//                probability +=
//                        hmm.getPi(i) *
//                                hmm.getOpdf(i).probability(oseq.get(0)) * beta[0][i];
        else
            for (int i = 0; i < hmm.nbStates(); i++) {
                probability = probability.add(BigDecimal.valueOf(hmm.getPi(i) * hmm.getOpdf(i).probability(oseq.get(0))).multiply(beta[0][i], m_Precision), m_Precision);
            }
    }

    public BigDecimal probabilityByBigDecimal() {
        return probability;
    }


}
