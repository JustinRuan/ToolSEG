//Copyright 2017 Jun Ruan
//
//        Licensed under the Apache License, Version 2.0 (the "License");
//        you may not use this file except in compliance with the License.
//        You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//        Unless required by applicable law or agreed to in writing, software
//        distributed under the License is distributed on an "AS IS" BASIS,
//        WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//        See the License for the specific language governing permissions and
//        limitations under the License.
package edu.whut.info.util;

import java.io.Serializable;


/**
 * @author Justin
 * @version 1.0
 * @created 08-11-2014 10:19:41
 */
public class GeneRange implements Serializable {

    /**
     *
     */
    private static final long serialVersionUID = -9151595175558457911L;
    public short Chr_id;
    public int Start;
    public int End;//not include the position of "End"

    public GeneRange() {
        Chr_id = 0;
        Start = 0;
        End = 0;
    }

    public GeneRange(short id, int s, int e) {
        Chr_id = id;
        Start = s;
        End = e;
    }

    public int length() {
        return End - Start;
    }

    public String toString() {
        return String.format("[%02d:% 6d - % 6d]", Chr_id, Start, End);
    }

}