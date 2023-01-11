### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 89c6c3ed-3129-4fe4-82cc-ec65d3f054db
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add(url = "https://github.com/TakekazuKATO/TailRec.jl")
	Pkg.add(["BenchmarkTools", "DataFrames", "PlutoUI", "PlutoPlotly", "Random"])
	using BenchmarkTools, TailRec, DataFrames, PlutoUI, PlutoPlotly, Random
end

# ╔═╡ bf238982-5291-4210-b055-211e80387e21
md"""
# Longest Increasing Subsequence Problem
Given a list of real numbers, e.g. {1, 5, 3, 10, 22},  consider subsequences of this list in the same order that are strictly increasing.  I.e. l[i+1] > l[i] for all i.  The question is what is the longest such subsequence of any arbitrary list provided.
"""

# ╔═╡ 2c1d5d9e-87f7-4861-b718-26e21c85e939
md"""
## Approach 1: Dynamic Programming
Take the input list as a vector $A$ and consider a new list of length equal to $A$ called $DP$.  This list will contain the length of the longest increasing subsequence ending at that point in $A$ i.e. $DP[i] = \text{length of longest increasing subsequence ending at } A[i]$.  So $DP[1] = 1$ always by definition.  Consider how to determine $DP[2]$. 
We can check to see if $A[2]$ can be appended to the end of any previous list.  Such a list would then be increased in length by 1.  The longest such list will then be the longest list ending at $A[2]$.  So we can iterate through all the longest lists found so far that terminate prior to $A[2]$ and only keep the longest resulting list of this comparison.  While this method of iteration will conveniently compare all possible subsequences, the number of checks that have to be done in general is $O(n^2)$.

With this method a list of about 100,000 numbers can be processed in ~7 seconds but if we increase the list length by 10x we would expect a run time 100x longer which is starting to get inconvenient.
"""

# ╔═╡ 2495339d-248c-421e-ab13-45cad597616b
@tailrec function lis_dp(A, i = 2, DP = fill(1, length(A)))
	#after checking each element return the list of longest sequences ending at each value
	i > length(A) && return (DP, findmax(DP))
	
	#update the longest sequence terminating at i if it can be appended and results in a list longer than any found so far
	for j = i-1:-1:1
		newl = DP[j] + 1 # length of new list if appened at index j
		(newl > DP[i]) && # ignore if shorter than previous best
		(A[j] < A[i]) && # ignore if cannot be appended
		(DP[i] = newl)
	end
	lis_dp(A, i+1, DP) #find longest subsequence of next element
end

# ╔═╡ 2370817f-914c-4644-90db-396425bf42cf
samplelist = [ 2, 5, 3, 7, 11, 8, 10, 13, 6 ]

# ╔═╡ 3ef2e6bc-4c4a-4fbe-a445-5c8946819ac1
lis_dp(samplelist)

# ╔═╡ d1dd8a4d-f45f-4cc2-aedc-1aedf9875851
A = rand(1:1_000_000, 100_000)

# ╔═╡ 69552382-f204-451f-a686-56254375bdb5
lis_dp(A)

# ╔═╡ ca817b74-4cc1-4701-b844-6d356c9ecb29
md"""
## Approach 2: Sorted Search
If we iterate through the list one element at a time one can notice that the longest possible subsequence observed up to that point can have a maximum length of i where i is the current index.  Furthermore, we can only append an element if it is larger than the end of an active list.  So if we have two active lists of equal length, the one with the smaller end value is the only one worth keeping.  At any given moment, we have a candidate subsequence that is the longest observed so far, but any other candidate shorter than this list is only worth considering if its ending value is smaller than the longer list.  Using these observations, we can design a method that keeps track of the ending values of candidate lists in increasing order of both list length and ending value.  Let's call this list listends with the property $\text{listends}[i] > \text{listends}[i-1]$.  Also, each time we encounter a new value, a list length can at most be increased by one, so the index i also represents the length of each candidate list as well as its ending value.

Initially, $\text{listends}[1] = A[1]$ representing the candidate list of length 1 that is simply the first element of A.  Now consider what happens at $i = 2$.  We observe $A[2]$ and have the following cases:

Case 1: $A[2] > \text{listends}[1]$

- In this case, $A[2]$ can be appended to the first candidate.  However we still should maintain the original candidate as future values could still be appended to it.  So listends will be transformed by $\text{listends}[2] = A[2]$.  So now $\text{listends} = {A[1], A[2]}$ representing two candidate lists of length 1 and 2.

Case 2: $A[2] \leq \text{listends}[1]$

- In this case, $A[2]$ cannot be appended to the first candidate but is an equal or better list of length 1 due to its value.  So we should replace the list of length 1 with this value and listends will be transformed by $\text{listends}[1] = A[2]$.  So now $\text{listends} = {A[2]}$ representing one new candidate list of length 1.

Imagine now we have many candidate lists and are considering a new index i.  Since listends is sorted we can efficiently find the index j such that $\text{listends}[k] \geq A[i] \forall k \geq j$.  $A[i]$ can only be appended to lists at index j-1 and lower and since all candidates less than j-1 are shorter lists, it is only worth considering the new list at j-1.  Again we consider the cases on j.

Case 1: $j = 1$

- In this case $A[i]$ is smaller than any candidate list end and is our new candidate for best list of length 1.  listends will be transformed by $\text{listends}[1] = A[i]$

Case 2: $j > \text{length}(\text{listends})$

- In this case $A[i]$ is larger than any candidate list end and will create a new list of $\text{length}(\text{listends}) + 1$.  listends will be transformed by $\text{listends}[j] = A[i]$

Case 3: $1 < j \leq \text{length}(\text{listends})$

- In this case $A[i]$ can be appended to the list at index j - 1, but its ending value will be less than or equal to the list currently at index j.  That means that it can replace that list, so listends will be transformed by $\text{listends}[j] = A[i]$

After examining these 3 cases we can observe that they are actually all the same regarding what action is taken to update listends.  Either an existing value is replaced at index j or a new value is appended to the end of the list.  This approach is implemented below in the function `lis`.
"""

# ╔═╡ 4058eb00-d518-497a-b1f4-3f971c9f81f4
@tailrec function lis(A::Vector{Int64}, i = 2, numlists = 1, listends = [A[1]])
	(i > length(A)) && return (listends = numlists, maxlength = numlists) #return listends and the list length after exhausting every element of A

	#find the first list end that is larger than A[i], note that listends is always in sorted order so this can be found efficiently with binary search
	j = searchsortedfirst(listends, A[i])

	if j <= numlists
		#case 1, current element is smaller than 1 or more of the current listends and thus can extend and replace an existing list
		listends[j] = A[i]
	else
		#case 2, current element is largest among all end candidates of active lists and should extend the number of active lists by 1
		push!(listends, A[i])
		numlists += 1
	end

	#check the next candidate element
	lis(A, i+1, numlists, listends)
end

# ╔═╡ f38d345d-4c70-4784-8af1-89165a1bcab8
md"""
Since the most candidate lists we could ever have is the length of A.  We can make a slight improvement by preallocating listends to its maximum length.  That way instead of ever appending elements we can simply always insert the element at the sorted index.  We must be careful though to only perform the sort operation on the list cut off at the appropriate listcount.  This improvement is implemented below as `lis_efficient`
"""

# ╔═╡ 71c15722-df8d-4cee-9bce-ce41c352219d
@tailrec function lis_efficient(A::Vector{Int64}, i = 2, numlists = 1, listends = fill(A[1], length(A)))
	(i > length(A)) && return (listends = numlists, maxlength = numlists)

	#Find the first list end that is larger than or equal to A[i].  This means that A[i] is strictly greater than the list end at j - 1 and could be added to that list. Note that listends is always in sorted order so the search can be done efficiently with binary search implemented here by `searchsortedfirst`
	j = searchsortedfirst(view(listends, 1:numlists), A[i])
	
	#insert in position j the list element at position i
	listends[j] = A[i]

	#if j > numlists we have added a new list, otherwise we are replacing an existing one.  `numlists + (j > numlists)` calculates this by abusing the integer value of booleans. 
	lis_efficient(A, i+1, numlists + (j > numlists), listends)
end

# ╔═╡ a13062e1-b16e-4e11-936c-45b62787985e
@btime lis($A)

# ╔═╡ 031234ed-7c63-4d0b-ba45-b8936d1fc511
@btime lis_efficient($A)

# ╔═╡ 483075a3-0d79-46f9-a9f7-1308a501186a
A_long = rand(1:100_000_000, 100_000_000)

# ╔═╡ adef84c0-09f8-4eb8-bf2c-51f81abee87b
lis(A_long)

# ╔═╡ d7c5ccd6-7ecf-4a79-8a74-b96d133dd38b
lis_efficient(A_long)

# ╔═╡ 8b78ce92-e9d6-4b27-83d3-8493d8511cdd
md"""
Note that these functions find the answer for the previous list A in about 3 ms vs 7 seconds for the dynamic programming method.  In fact, even for a list of length 100 million, they can run in ~5 seconds.  The advantage of these functions is that for each of the n elements of A, they need to perform a search operation on a sorted list of maximum length i where i < n.  Since the list is sorted, such a search can be done with binary search reducing the number of operations to the order of log(i).  It isn't quite as good as linear scaling, but it makes dealing with extremely long lists tractable.
"""

# ╔═╡ 2443519c-d39e-40c7-b0ac-0ffe4b1c950f
md"""
# Modification 1: Sum of Longest Increasing Subsequence
"""

# ╔═╡ 128502dc-2cc5-40fe-bab5-e6c3cec04ba8
md"""
Consider the same problem with one modification.  Instead of finding the subsequence with the longest length, we are interested in the one with the largest sum.  

## Approach 1: Dynamic Programming
Due to the nature of the original dynamic programming solution iterating through all subsequences, it can be easily modified to suit this new objective.  Our array now stores sum of the subsequence ending at each element with the highest value rather than the list length and everything else is unchanged.  This algorithm is implemented as `lis_sum_dp` below.
"""

# ╔═╡ 8e73e8de-4906-4d6d-b92e-f3bd6684b57a
@tailrec function lis_sum_dp(A, i = 2, DP = copy(A))
	#after checking each element return the list of highest sum sequences ending at each value
	i > length(A) && return (DP, findmax(DP))
	
	#update the longest sequence terminating at i if it can be appended and results in a list longer than any found so far
	for j = i-1:-1:1
		newsum = DP[j] + A[i] # sum of new list if appened at index j
		(newsum > DP[i]) && # ignore if smaller than previous best
		(A[j] < A[i]) && # ignore if cannot be appended
		(DP[i] = newsum)
	end
	lis_sum_dp(A, i+1, DP) #find longest subsequence of next element
end

# ╔═╡ 2e0b40a5-04b1-4b87-bfa1-8abc4ef244ed
B = [100, 3, 4, 5, 10]

# ╔═╡ 99440ab2-3f84-4232-b712-66da4e6fbb1f
lis_sum_dp(B)

# ╔═╡ 2ecda151-3bec-47e2-8763-6f4a14b9cf49
md"""
This method has the same computational complexity as the previous DP solution which makes it intractable for very long lists.  The question remains again, can we apply a more efficient method to solve this problem similar to the solution for list length?  Luckily it is possible to use some of the previous techniques on this problem but the modifications are not as trivial as the DP case.
"""

# ╔═╡ 58ed6bc4-48ac-459f-b50a-e1003cafd1d7
md"""
## Approach 2: Sorted Search
"""

# ╔═╡ 0cece5b1-7762-4ad2-a45e-aa7027820fe5
#how can this be adapted to instead find the strictly increasing sequence with the largest sum?  In this solution we keep active lists such that the end of each list of increasing length increases.  This is consistent with the criteria that a shorter sequence only has the potential to eventually become longer if its end value is smaller than any sequence of larger size.  In the past if we find a sequence of length 1 whose value is smaller than the existing list of length 1 then it should replace it, however, for the purpose of the sum, it might still be useful to keep the larger element.  So it seems we must keep both.  What if we maintain an array S[j] where S[j] is the largest sum so far for any subsequence ending at index j?  If we have two lists that share the same sum, we should only keep the only that ends in the smaller element.  So we need to maintain an array of sums and corresponding ending elements.  We should have two sorting orders for it, one in terms of sum and the other in terms of ending element.  When we look at a new element we can only add it to the lists that are above the midpoint value and then these will have a new ending value.  Then we need to make sure that we don't have any redundant lists

# ╔═╡ 536b73e1-3d5a-4269-8504-a2f19a88f0d2
#Consider the following example A = [ 2, 5, 3, 7, 11, 8, 10, 13, 6 ]
#Starting with the first element we have a sequence ending in 2 with a sum of 2 so (2, 2)
# Now we encounter a 5.  It could be added so now we have (2, 2), (5, 7).  Appending it like this maintains the order of both sum and ending value.

# Now we encounter a 3.  It can only be appended to the first list so we have:
# (2, 2), (3, 5), (5, 7).  Seems like we should keep all 3.

# Now we encounter a 7.  It can be appened to all 3 lists so this will double the number of active lists: (2, 2), (3, 5), (5, 7), (7, 9), (7, 12), (7, 14).  So we appened 7 to every list, however, it is only worth keeping the one of these with the highest sum so we actually only need to append 7 to the end of the ending pair resulting in: (2, 2), (3, 5), (5, 7), (7, 14)

# After adding 11 we have a similar situation again: (2, 2), (3, 5), (5, 7), (7, 14), (11, 25)

# The next element is 8, by the same logic as before, this should be added to (5, 7) and inserted in there yielding: (2, 2), (3, 5), (5, 7), (8, 22), (11, 25)

# The next element is 10, this should be added to (8, 22) yielding: (2, 2), (3, 5), (5, 7), (8, 22), (10, 32), (11, 25).

#Now for the first time we have a situation where there is a list that ends in 10 but has a higher sum than a list that ends in 11.  This means that we should eliminate the list that currently is at the end.  In general if we insert a new list and the sum value ends up exceeding any list above it those lists should be eliminated.

# ╔═╡ 036e49cb-6b6e-4bdf-9975-42c58a74ebba
#this solution is better than the dynamic programming version but ideally it would be done without the insert! and delete! operations and just having an array preallocated of the maximum length and just shifting data around in it. 
@tailrec function lis_sum(A::Vector{Int64}, i = 2, listends = [A[1]], listsums = [A[1]])
	(i > length(A)) && return listsums[end]

	numlists = length(listends)

	#Find the first list end that is larger than or equal to A[i].  This means that A[i] is strictly greater than the list end at j - 1 and could be added to that list.  If j is 1 then A[i] is smaller than any list end and thus cannot be added to any list.  If j is numlists + 1 then A[i] is larger than any list end and thus can be added to the longest list so far.  If j is in between then there is an existing list of length j that should be replaced by adding A[i] to the previous list ending in A[j-1].  Note that listends is always in sorted order so the search can be done efficiently with binary search implemented here by `searchsortedfirst`
	j = searchsortedfirst(listends, A[i])
	
	#calculate the sum of the new list that will be at position j.  If j = 1 then this is a new list with no previous elements
	newsum = A[i] + (j > 1 && listsums[j-1])

	#check to see if any items should be discarded from the existing lists.  Prior to inserting only elements from k + j - 1 onward should be kept.  In the case where we keep everything, k = 1.  Thus the number of eliminated values is k-1
	k = searchsortedfirst(view(listsums, j:numlists), newsum, lt = <=)

	#delete values from listsums and listends that are redundant note that if k = 1, nothing happens, if k = 2, then 1 element is deleted etc...
	for ind in 1:k-1
		deleteat!(listsums, j)
		deleteat!(listends, j)
	end
	
	#insert in position j the list element at position i
	insert!(listends, j, A[i])
	
	#insert the updated sum here as well
	insert!(listsums, j, newsum)
	
	#if j > 1 we have added a new list, otherwise we are replacing an existing one
	lis_sum(A, i+1, listends, listsums)
end

# ╔═╡ efa4e643-d495-486c-a33c-b04a936ad89f
#this solution is better than the dynamic programming version but ideally it would be done without the insert! and delete! operations and just having an array preallocated of the maximum length and just shifting data around in it. 
@tailrec function lis_sum_efficient(A::Vector{Int64}, i = 2, numlists = 1, listends = fill(A[1], length(A)), listsums = fill(A[1], length(A)))
	(i > length(A)) && return listsums[numlists]

	#Find the first list end that is larger than or equal to A[i].  This means that A[i] is strictly greater than the list end at j - 1 and could be added to that list.  If j is 1 then A[i] is smaller than any list end and thus cannot be added to any list.  If j is numlists + 1 then A[i] is larger than any list end and thus can be added to the longest list so far.  If j is in between then there is an existing list of length j that should be replaced by adding A[i] to the previous list ending in A[j-1].  Note that listends is always in sorted order so the search can be done efficiently with binary search implemented here by `searchsortedfirst`
	j = searchsortedfirst(view(listends, 1:numlists), A[i])
	
	#calculate the sum of the new list that will be at position j.  If j = 1 then this is a new list with no previous elements
	newsum = A[i] + (j > 1 && listsums[j-1])

	#check to see if any items should be discarded from the existing lists.  Prior to inserting only elements from k + j - 1 onward should be kept.  In the case where we keep everything, k = 1.  Thus the number of eliminated values is k-1 and the number of kept values is numlists - j + 2 - k
	k = searchsortedfirst(view(listsums, j:numlists), newsum, lt = <=)


	listsums[j+1:numlists+2-k] .= view(listsums, j-1+k:numlists)
	listends[j+1:numlists+2-k] .= view(listends, j-1+k:numlists)
	
	#insert in position j the list element at position i
	listends[j] = A[i]
	
	#insert the updated sum here as well
	listsums[j] = newsum
	
	#if j > 1 we have added a new list, otherwise we are replacing an existing one
	lis_sum_efficient(A, i+1, numlists + 2 - k, listends, listsums)
end

# ╔═╡ db268908-e2c7-4666-a7a1-e46800b00d20
lis_sum(B)

# ╔═╡ 776d70b6-504f-4e00-9c1a-f00a336e3da4
lis_sum_dp(A)

# ╔═╡ ca934474-ca2e-4ecf-ba3e-8178f436154a
lis_sum(A)

# ╔═╡ 66e85ede-1962-4d5c-99fe-25478a0536af
lis_sum_efficient(A)

# ╔═╡ 21576f52-95f6-4a63-b3fe-07b3944b1211
md"""
# 2D Extension: Tallest Box Stack Problem
"""

# ╔═╡ df5304f7-7b50-40d3-b38a-75ddac8e139d
getdim(n = 1_000_000) = rand(1:n)

# ╔═╡ 5fd9bbde-8cc4-11ed-084c-938430ee369a
function makebox(n = 1_000_000)
	keys = [:l, :w, :h]
	values = [getdim(n) for _ in 1:3]
	(; zip(keys, values)...)
end

# ╔═╡ c5d0eb69-cc7b-4c24-8337-579e58bb1e66
function makeboxes(n, d = n)
	df = DataFrame(makebox(d) for _ in 1:n)
	sort(select(df, :l, :w, :h, [:l, :w] => ByRow((a, b) -> a * b) => :area), :area)
end

# ╔═╡ f165186d-659d-4792-bddc-6f2f3d526098
#checks if b2 can be stacked on b1 or if b1 can be stacked under b
function canstack(b1, b2)
	(b2[1] < b1[1]) && (b2[2] < b1[2])
end

# ╔═╡ c4f4ec52-ed81-47f3-a1cb-e3bc42860e8f
boxes = makeboxes(10)

# ╔═╡ 331f36f2-860b-4a94-a442-70b7c255c391
boxes.h

# ╔═╡ 9862b027-989b-4f7c-aa83-65452928888a
canstack(boxes[1, :], boxes[2, :])

# ╔═╡ fc895b19-b577-442c-b14d-33fec233d18a
function tallest_stack_simple(boxes)
#iterate through every posssible stack in n^2 complexity assuming boxes are ordered by area
	numboxes = size(boxes, 1)
	heights = copy(boxes.h) #should be the highest stack with a base at this index
	stackinds = [[i] for i in 1:numboxes]

	for i = 2:numboxes
   		for j = i - 1:-1:1	
		  	if canstack(boxes[i, :], boxes[j, :]) 
			  	newheight = boxes[i, :h] + heights[j]
			  	if (newheight > heights[i]) 
					heights[i] = newheight
					stackinds[i] = vcat(stackinds[j], i)
				end
			end
		end
	end
	return (stackheights = heights, stacks = stackinds, bestheight = maximum(heights))
end

# ╔═╡ 7021121c-0c6d-4d0e-96ac-fe65b83d4e6a
#Now let's apply this to the stacked box problem, assume that boxes is a dataframe of boxes containing length, width, height, and the area of each box.  Sorted by area from smallest first.  Note that a box can only successfully placed below another if it has strictly larger area.  Even if this is the case additional checks on length and width must also be made.
@tailrec function tallest_stack(boxes, i = 2, bases = [boxes[1, :]], areas = [boxes[1, :].area], heights = [boxes[1, :].h])
	(i > size(boxes, 1)) && return (base_areas = areas, stackheights = heights, bestheight = maximum(heights))

	checkbox(base, newbox) = (newbox.area > base.area) && (newbox.w > base.w) && (newbox.l > base.l)

	numlists = length(heights)

	#We maintain two lists of equal length.  bases contains the current base box of all active stacks.  Heights contains the current height of every active stack.  At the beginning we have a single element for each list consisting of the first box in the list which is sorted by base area.  If we add a new box to active stacks then those new stacks will all have the same larger base area and appear identical with respect to adding additional boxes.  Therefore, only the tallest of such stacks needs to be kept at any given moment.  Similarly consider the base of any two stacks in which one base is of larger area than another.  That base will be strictly more difficult to add to than the other one, therefore, such a stack is only worth considering if it is taller.  Using this method of addition we maintain both lists such that every active stack has increasing base areas and heights in tandem.  Each new box added will create potentially one more active stack and potentially remove other active stacks if they have larger base areas but lower heights.

	#We could encounter a new box with an equal base area to the previous base in which case we cannot stack it but it could be potentially stacked on any of the smaller bases to a point yielding a better stack.  

	#check in a range of indices going backwards the first base that is stackable.  In the worst case scenario we will have to check every single previous active stack.  Is there any way to speed up this step?  We know that this newbox has an area greater than or equal to all previous stacks.  In the existing list we can already exclude any base with an equal area so one idea is to also have a sorted list in terms of area and use that to shorten the amount of checks.
	@tailrec function checkbases(j, minindex)
		j < minindex && return 0
		checkbox(bases[j], boxes[i, :]) && return j
		return checkbases(j-1, minindex)
	end

	#w1*l1 > w2*l2 && w1 > w2 what does this imply about l1 and l2?  l1/l2 > w2/w1
	#For example let's say we have something that's 5x5 and another thing that is 2x30 so we have both 60>25 and 30>5 yet it is not stackable.  If we maintain two sorted lists in terms of length and width then we would have 2 sets of indices of potential stacks and the most we'd have to check of each is the minimum of the two sets of indices since we take the intersection.  If this number is small enough then it would be worth iterating through that and picking the one with the tallest stack rather than going through each one.

	function updatestacks!(newbox)
		j = checkbases(numlists, 1)

		#if j is 0 then the new stack height will simply be the current box
		topheight = j == 0 ? 0 : heights[j]

		#calculate the height of the new stack
		newheight = newbox.h + topheight

		#find the appropriate slot for the new stack
		k = searchsortedfirst(heights, newheight)
	
		insert!(heights, k, newheight)
		insert!(areas, k, newbox.area)
		insert!(bases, k, newbox)
	end
	updatestacks!(boxes[i, :])
	
	tallest_stack(boxes, i+1, bases, areas, heights)
end

# ╔═╡ afc88511-6770-4dc9-9b5d-e3a5d6472c88
testboxes = makeboxes(100_000)

# ╔═╡ 98d743df-83e5-4a61-8ab5-b04c30b17e8c
tallest_stack(testboxes)

# ╔═╡ 81f370e8-7cff-4da1-acbf-0a29d738d4e9
# ╠═╡ disabled = true
#=╠═╡
tallest_stack_simple(testboxes)
  ╠═╡ =#

# ╔═╡ 0328d96c-106d-456e-8023-e6dae14a37f5
#I suspect this doesn't actually work because the lookup gets screwed up every time something new is added.
#Now let's apply this to the stacked box problem, assume that boxes is a dataframe of boxes containing length, width, height, and the area of each box.  Sorted by area from smallest first.  Note that a box can only successfully placed below another if it has strictly larger area.  Even if this is the case additional checks on length and width must also be made.
@tailrec function tallest_stack_efficient(boxes, thresh = 0.1, i = 2, heightinds = [1], heightindslookup = [1], areas = [boxes[1, :].area], heights = [boxes[1, :].h], lengths = [boxes[1, :l]], lengthinds = [1], widths = [boxes[1, :w]], widthinds = [1])
	
	(i > size(boxes, 1)) && return (base_areas = areas, stackheights = heights, bestheight = maximum(heights))

	checkbox(base, newbox) = (newbox.area > base.area) && (newbox.w > base.w) && (newbox.l > base.l)

	numlists = length(heights)

	l_check = searchsortedfirst(lengths, boxes[i, :l])
	w_check = searchsortedfirst(widths, boxes[i, :w])

	maxvalid = min(l_check - 1, w_check - 1)

	@tailrec function checkbases2(checkinds, j = 1, bestheight=boxes[i, :h], bestind = 0)
		j > length(checkinds) && return (bestheight, bestind)
		if checkbox(boxes[heightindslookup[checkinds[j]], :], boxes[i, :])
			newheight = heights[heightindslookup[checkinds[j]]] + boxes[i, :h]
			if newheight > bestheight
				bestheight = newheight
				bestind = j
			end
		end
		return checkbases2(checkinds, j+1, bestheight, bestind)
	end	

	#look through the bases starting at the maximum height and going backwards to find the first one that stacks
	@tailrec function checkbases(j, minindex)
		j < minindex && return 0
		checkbox(boxes[heightinds[j], :], boxes[i, :]) && return j
		return checkbases(j-1, minindex)
	end

	function updatestacks!(newbox)
		j = checkbases(numlists, 1)

		#if j is 0 then the new stack height will simply be the current box
		topheight = j == 0 ? 0 : heights[j]

		#calculate the height of the new stack
		newheight = newbox.h + topheight

		#find the appropriate slot for the new stack
		k = searchsortedfirst(heights, newheight)
	
		insert!(heights, k, newheight)
		insert!(areas, k, newbox.area)
		insert!(heightinds, k, i)
		insert!(heightindslookup, i, k)
	end

	if maxvalid < numlists * thresh
		checkinds = intersect(lengthinds[1:l_check-1], widthinds[1:w_check-1])
		(bestheight, bestind) = checkbases2(checkinds)
		#find the appropriate slot for the new stack
		k = searchsortedfirst(heights, bestheight)
		insert!(heights, k, bestheight)
		insert!(areas, k, boxes[i, :area])
		insert!(heightinds, k, i)
		insert!(heightindslookup, i, k)
	else
		updatestacks!(boxes[i, :])
	end
	
	insert!(lengths, l_check, boxes[i, :l])
	insert!(lengthinds, l_check, i)
	insert!(widths, w_check, boxes[i, :w])
	insert!(widthinds, w_check, i)
	
	
	tallest_stack_efficient(boxes, thresh, i+1, heightinds, heightindslookup, areas, heights, lengths, lengthinds, widths, widthinds)
end

# ╔═╡ dc49899c-ce32-4d79-9a69-82e8c00db45f
testboxes2 = makeboxes(100_000)

# ╔═╡ 53419468-9793-4186-a8b2-263dcbf606df
tallest_stack(testboxes2)

# ╔═╡ 9c67be56-79bf-450d-a094-5b5b4e1d06b6
tallest_stack_efficient(testboxes2, 0.25)

# ╔═╡ 3a67e330-845d-4b48-93b1-d47e3c0ae026
PlutoUI.TableOfContents()

# ╔═╡ Cell order:
# ╠═89c6c3ed-3129-4fe4-82cc-ec65d3f054db
# ╟─bf238982-5291-4210-b055-211e80387e21
# ╟─2c1d5d9e-87f7-4861-b718-26e21c85e939
# ╠═2495339d-248c-421e-ab13-45cad597616b
# ╠═2370817f-914c-4644-90db-396425bf42cf
# ╠═3ef2e6bc-4c4a-4fbe-a445-5c8946819ac1
# ╠═d1dd8a4d-f45f-4cc2-aedc-1aedf9875851
# ╠═69552382-f204-451f-a686-56254375bdb5
# ╟─ca817b74-4cc1-4701-b844-6d356c9ecb29
# ╠═4058eb00-d518-497a-b1f4-3f971c9f81f4
# ╟─f38d345d-4c70-4784-8af1-89165a1bcab8
# ╠═71c15722-df8d-4cee-9bce-ce41c352219d
# ╠═a13062e1-b16e-4e11-936c-45b62787985e
# ╠═031234ed-7c63-4d0b-ba45-b8936d1fc511
# ╠═483075a3-0d79-46f9-a9f7-1308a501186a
# ╠═adef84c0-09f8-4eb8-bf2c-51f81abee87b
# ╠═d7c5ccd6-7ecf-4a79-8a74-b96d133dd38b
# ╟─8b78ce92-e9d6-4b27-83d3-8493d8511cdd
# ╟─2443519c-d39e-40c7-b0ac-0ffe4b1c950f
# ╟─128502dc-2cc5-40fe-bab5-e6c3cec04ba8
# ╠═8e73e8de-4906-4d6d-b92e-f3bd6684b57a
# ╠═2e0b40a5-04b1-4b87-bfa1-8abc4ef244ed
# ╠═99440ab2-3f84-4232-b712-66da4e6fbb1f
# ╟─2ecda151-3bec-47e2-8763-6f4a14b9cf49
# ╟─58ed6bc4-48ac-459f-b50a-e1003cafd1d7
# ╠═0cece5b1-7762-4ad2-a45e-aa7027820fe5
# ╠═536b73e1-3d5a-4269-8504-a2f19a88f0d2
# ╠═036e49cb-6b6e-4bdf-9975-42c58a74ebba
# ╠═efa4e643-d495-486c-a33c-b04a936ad89f
# ╠═db268908-e2c7-4666-a7a1-e46800b00d20
# ╠═776d70b6-504f-4e00-9c1a-f00a336e3da4
# ╠═ca934474-ca2e-4ecf-ba3e-8178f436154a
# ╠═66e85ede-1962-4d5c-99fe-25478a0536af
# ╟─21576f52-95f6-4a63-b3fe-07b3944b1211
# ╠═df5304f7-7b50-40d3-b38a-75ddac8e139d
# ╠═5fd9bbde-8cc4-11ed-084c-938430ee369a
# ╠═c5d0eb69-cc7b-4c24-8337-579e58bb1e66
# ╠═f165186d-659d-4792-bddc-6f2f3d526098
# ╠═c4f4ec52-ed81-47f3-a1cb-e3bc42860e8f
# ╠═331f36f2-860b-4a94-a442-70b7c255c391
# ╠═9862b027-989b-4f7c-aa83-65452928888a
# ╠═fc895b19-b577-442c-b14d-33fec233d18a
# ╠═7021121c-0c6d-4d0e-96ac-fe65b83d4e6a
# ╠═afc88511-6770-4dc9-9b5d-e3a5d6472c88
# ╠═98d743df-83e5-4a61-8ab5-b04c30b17e8c
# ╠═81f370e8-7cff-4da1-acbf-0a29d738d4e9
# ╠═0328d96c-106d-456e-8023-e6dae14a37f5
# ╠═dc49899c-ce32-4d79-9a69-82e8c00db45f
# ╠═53419468-9793-4186-a8b2-263dcbf606df
# ╠═9c67be56-79bf-450d-a094-5b5b4e1d06b6
# ╠═3a67e330-845d-4b48-93b1-d47e3c0ae026
