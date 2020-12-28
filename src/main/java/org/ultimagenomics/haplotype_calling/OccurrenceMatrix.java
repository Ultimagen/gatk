package org.ultimagenomics.haplotype_calling;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import java.util.*;
import java.util.stream.Collectors;

public class OccurrenceMatrix<R,C> {

    private List<R> rowNames;
    private List<C> colNames;
    private Map<C, Integer> col2idx;
    private int nRows;
    private int nCols;
    private boolean[][] occurrenceMatrix;

    public OccurrenceMatrix(final Map<R, Collection<C>> input_map){
        nRows = input_map.size();
        rowNames = input_map.keySet().stream().collect(Collectors.toList());

        col2idx = new HashMap<>();
        int col_count = 0;
        for ( R row: input_map.keySet() ){
            for (C col: input_map.get(row)){
                if (!col2idx.containsKey(col)){
                    col2idx.put(col, col_count);
                    col_count++;
                }
            }
        }
        nCols = col_count;
        colNames = new ArrayList<>();
        for ( C col: col2idx.keySet() ){
            colNames.add(col2idx.get(col), col);
        }


        occurrenceMatrix = new boolean[nRows][nCols];
        for (int r = 0; r < nRows; r++){
            for ( C col: input_map.get(rowNames.get(r)))
                occurrenceMatrix[r][col2idx.get(col)] = true;
        }
    }

    public List<Pair<C, C>> nonCoOcurringColumns() {

        List<Pair<Integer, Integer>> result = new ArrayList<>();
        for (int i = 0; i < nCols; i++) {
            for (int j = i + 1; j < nCols; j++) {
                boolean flag = false;
                for (int r = 0; r < nRows; r++) {
                    if (occurrenceMatrix[r][i] & occurrenceMatrix[r][j]) {
                        flag = true;
                        break;
                    }
                }
                if (!flag) {
                    result.add(ImmutablePair.of(i,j));
                }

            }
        }
        List<Pair<C,C>> vc_result = new ArrayList<>();
        for (Pair<Integer, Integer> res: result) {
            vc_result.add(ImmutablePair.of(colNames.get(res.getLeft()), colNames.get(res.getRight())));
        }
        return vc_result;
    }

    public List<List<C>> getIndependentSets(){
        List<Pair<C,C>> nonCoOcurring = nonCoOcurringColumns();
        return null;
    }
}