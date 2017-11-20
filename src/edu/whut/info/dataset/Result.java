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
package edu.whut.info.dataset;

/**
 * Created by Justin on 2016/10/25.
 */
public class Result {
    public int id;
    public int pos; //location of breakpoint
    public int[] range;// continuous area
    public boolean isBreakPoint;
    public double value1;
}
