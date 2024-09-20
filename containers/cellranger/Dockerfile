FROM ubuntu:22.04
RUN apt-get update
RUN apt-get install -y wget
RUN wget -O cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1726860797&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=Jk8r2tpwMK0YvRei8i3m4Oh6xFjDXmk~qM3n9Nb9rFpWLoyMyvKrnY4RybiMj45xJdBoBgaQIPF7baJaJguZ78mmDdJ~64Pn7LdHkqi3cOhMnANTuBsNJgJJq1ryPSc2xUf31RtWZYgrTxcGuGrm~aljJ2MQxqL6nQxnjsp5FhiJblkXTVZBzk7qnhY5yTNWOFcqmighFH54O4CvTrCWbyxR3smwgMe3XnbBy2-Z2AIJzESX7DViivVEkriX0UaydAFB5u9s1wPqBYP~i4bQL4cwAVqKHb69DhxNFpx7Sp-BqRIQL~-SraS3Y5zOe-lsKdxS~nY2dGskEBrdgRbkKQ__"
RUN tar -xvf cellranger-8.0.1.tar.gz
ENV PATH=./cellranger-8.0.1:$PATH
